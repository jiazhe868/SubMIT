#!/bin/bash
# step3_do_inversions.sh - run the MCMC inversions (finv) for every number of subevents.
#
# Without arguments: every inv_<event>_<N>sub folder, N = 1, 2, ... in order. After N = 1 the
# 1-subevent calibration (prepare_SubMIT/calibrate_from_1sub.sh) tightens the bounds and windows of
# the N >= 2 runs, balances the wave-type weights and screens stations against the 1-subevent
# prediction; N = 1 is then run again under that same configuration, so every point of the L-curve
# is scored on identical data (otherwise N = 1 is scored on longer windows and other weights and the
# 5% rule compares incomparable misfits). The first N = 1 result is kept as calibration_pass_model.txt
# (a name that the *best.dat globs below cannot match: as calibration_pass_best.dat it was deleted by
# the rerun's clean-up).
# With arguments: only those N, without calibration (grid extensions, reruns after editing).
#
# Warm start: before N >= 2 runs, the best model of N-1 is copied to warm_start.dat; half of the
# chains start from it plus one subevent drawn from the prior (finv.f90).
#
# History: this replaces step3_do_inversions.sh + step3_do_inversions_tempnsub.sh, which ran each
# folder through awk system() (killing the driver orphaned the awk, which kept starting mpirun jobs).
#
# usage (in IRIS/):  bash ../programs/step3_do_inversions.sh          all N, with calibration
#                    bash ../programs/step3_do_inversions.sh 6 7      only N = 6 and 7
# stop a run:        sh ../programs/stop_step3.sh "$PWD"
# ranks: SUBMIT_NCORE (default physical cores, at most 32); GNU toolchain: SUBMIT_NO_INTEL=1
if [ -z "$BASH_VERSION" ]; then exec bash "$0" "$@"; fi
# The whole script is one { ...; exit; } block, so bash parses all of it before running any of it.
# bash otherwise reads a script while it runs it, and an edit during a multi-hour step 3 (as the
# README invites for Par.file and the scripts) made it resume at a stale byte offset in the edited
# file: "syntax error: unexpected end of file" after N = 5 of a test run, before the L-curve.
{
set -u
# capture the arguments first: sourcing Intel's setvars.sh rewrites the positional parameters
requested_counts="$*"
set --
. ../programs/compiler_environment.sh

# Every station distance must have Green's functions: a missing file was once read as silent garbage
# for one station (2024 Chile AU.MAW). Check before spending hours of compute.
check_green_function_coverage() {
    local event=$1 inversion_dir="inv_${1}_1sub" gf_folder first_depth_folder station_list distance missing=0
    gf_folder="gf_${event}/greenFuncDir_disp"
    [ -d "$inversion_dir" ] && [ -d "$gf_folder" ] || return 0
    first_depth_folder=$(ls -d "$gf_folder"/vmodel_* 2>/dev/null | head -1)
    [ -n "$first_depth_folder" ] || return 0
    for station_list in stations.info stationsSH.info; do
        [ -s "$inversion_dir/$station_list" ] || continue
        for distance in $(gawk '{print $NF}' "$inversion_dir/$station_list" | sort -u); do
            if ! ls "$first_depth_folder/$distance".grn.* > /dev/null 2>&1; then
                echo "ERROR: no Green's functions at distance $distance (listed in $inversion_dir/$station_list)"
                missing=1
            fi
        done
    done
    return $missing
}

# build and run one inversion folder; best_model_hybrid.dat = best sample over all chains.
# mpirun reads /dev/null: it forwards its stdin to rank 0, and the stdin of this loop is
# event_list.dat, whose remaining events it would otherwise consume.
run_inversion() {
    local inversion_dir=$1
    ( cd "$inversion_dir" && make clean > /dev/null 2>&1 && make finv > make.log 2>&1 ) \
        || { echo "COMPILE FAILED in $inversion_dir (see $inversion_dir/make.log)"; return 1; }
    ( cd "$inversion_dir" && rm -f ./*best.dat best_model_*.dat ./*chain.dat mpi.out island_best_misfits.txt \
        && mpirun -np "$submit_rank_count" ./finv > mpi.out 2>&1 < /dev/null ) \
        || { echo "INVERSION FAILED in $inversion_dir (mpirun exit status; see $inversion_dir/mpi.out)"; return 1; }
    if ! ls "$inversion_dir"/[0-9]*best.dat > /dev/null 2>&1; then
        echo "INVERSION PRODUCED NO OUTPUT in $inversion_dir (see $inversion_dir/mpi.out)"
        return 1
    fi
    sort -gk3 "$inversion_dir"/[0-9]*best.dat | head -1 > "$inversion_dir/best_model_hybrid.dat"
    echo "inversion done: $inversion_dir best misfit $(gawk '{print $3}' "$inversion_dir/best_model_hybrid.dat")"
}

while read -r event_name rest; do
    check_green_function_coverage "$event_name" || { echo "step3 ABORTED: Green's-function coverage incomplete"; exit 1; }
    if [ -n "$requested_counts" ]; then
        subevent_counts=$requested_counts
    else
        subevent_counts=$(ls -d inv_"${event_name}"_*sub 2>/dev/null | sed 's/.*_\([0-9]*\)sub$/\1/' | sort -n)
    fi
    for subevent_count in $subevent_counts; do
        inversion_dir="inv_${event_name}_${subevent_count}sub"
        [ -d "$inversion_dir" ] || { echo "skipping $inversion_dir (missing)"; continue; }
        # warm start: half of the chains begin from the best N-1 model plus one subevent drawn from the
        # prior (finv.f90); each N refines the previous one instead of searching from scratch
        rm -f "$inversion_dir/warm_start.dat"
        previous_best="inv_${event_name}_$((subevent_count - 1))sub/best_model_hybrid.dat"
        [ "$subevent_count" -gt 1 ] && [ -s "$previous_best" ] && cp "$previous_best" "$inversion_dir/warm_start.dat"
        # a failed N stops step 3: the L-curve would otherwise select among the N that happened to work
        run_inversion "$inversion_dir" || { echo "step3 ABORTED at $inversion_dir"; exit 1; }
        if [ -z "$requested_counts" ] && [ "$subevent_count" -eq 1 ] && [ -f "inv_${event_name}_2sub/Par.file" ]; then
            bash ../programs/prepare_SubMIT/calibrate_from_1sub.sh "$inversion_dir" "$event_name" \
                || { echo "step3 ABORTED: the 1-subevent calibration failed (messages above)"; exit 1; }
            cp "$inversion_dir/best_model_hybrid.dat" "$inversion_dir/calibration_pass_model.txt"
            cp "inv_${event_name}_2sub/Par.file" "$inversion_dir/Par.file"
            echo "re-running $inversion_dir under the calibrated N>=2 configuration"
            run_inversion "$inversion_dir" || { echo "step3 ABORTED at the rerun of $inversion_dir"; exit 1; }
        fi
    done
done < event_list.dat
exit 0
}
