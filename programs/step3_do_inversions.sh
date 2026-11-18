#!/bin/bash
# step3_do_inversions.sh - run the MCMC inversions (finv) for every number of subevents.
#
# Without arguments: every inv_<event>_<N>sub folder, N = 1, 2, ... in order. After N = 1 the
# 1-subevent calibration (prepare_SubMIT/calibrate_from_1sub.sh) tightens the bounds and windows of
# the N >= 2 runs, balances the wave-type weights and screens stations against the 1-subevent
# prediction; N = 1 is then run again under that same configuration, so every point of the L-curve
# is scored on identical data (otherwise N = 1 is scored on longer windows and other weights and the
# selection of plot_lcurve.py compares incomparable misfits). The first N = 1 result is kept as
# calibration_pass_model.txt (a name that the *best.dat globs below cannot match: as
# calibration_pass_best.dat it was deleted by the rerun's clean-up).
# With arguments: only those N, without calibration (grid extensions, reruns after editing).
#
# After each inversion, prepare_SubMIT/evaluate_misfit_components.sh splits the best misfit of every
# chain into the data misfit and the two penalty factors (misfit_components.txt, drawn by
# plot_lcurve.py next to the L-curve; about 20 s per folder).
#
# Centroid range: the range of centroid times of every N comes from the magnitude (duration + 6 s,
# get_search_par.py). When the first 1-subevent run ends on its upper bound, the moment release
# lies later than that: the range of every N is doubled, the windows are lengthened by the same
# amount, and N = 1 runs again, at most twice (widen_centroid_range_while_on_bound below).
#
# Warm start: before N >= 2 runs, the best model of N-1 is copied to warm_start.dat; half of the
# chains start from it plus one subevent, drawn from the prior or, when that model leaves no free
# slot in time, inserted between two of its subevents (finv.f90), and the line after "inversion
# done" says how many did.
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
    ( cd "$inversion_dir" && rm -f ./*best.dat best_model_*.dat ./*chain.dat mpi.out island_best_misfits.txt misfit_components* \
        && mpirun -np "$submit_rank_count" ./finv > mpi.out 2>&1 < /dev/null ) \
        || { echo "INVERSION FAILED in $inversion_dir (mpirun exit status; see $inversion_dir/mpi.out)"; return 1; }
    if ! ls "$inversion_dir"/[0-9]*best.dat > /dev/null 2>&1; then
        echo "INVERSION PRODUCED NO OUTPUT in $inversion_dir (see $inversion_dir/mpi.out)"
        return 1
    fi
    sort -gk3 "$inversion_dir"/[0-9]*best.dat | head -1 > "$inversion_dir/best_model_hybrid.dat"
    echo "inversion done: $inversion_dir best misfit $(gawk '{print $3}' "$inversion_dir/best_model_hybrid.dat")"
    # none when the N-1 model leaves no room in time for another subevent (finv.f90): the search for
    # this N then started from the prior, and its misfit is less likely to improve on N-1
    if [ -s "$inversion_dir/warm_start.dat" ]; then
        echo "    $(grep -ac 'Like_init.* warm' "$inversion_dir/mpi.out") of $submit_rank_count chains started from the best model with one subevent fewer"
    fi
    # a diagnostic: it warns on failure but never stops step 3
    bash ../programs/prepare_SubMIT/evaluate_misfit_components.sh "$inversion_dir" \
        || echo "WARNING: misfit components failed for $inversion_dir"
    return 0
}

# A best 1-subevent model whose centroid time sits on the upper bound of its range was cut off by
# that bound. 2026 Colombia Mw 7.4, 116 km deep (v3, two seeds): every chain ended at 30.9-31.0 s of
# a range that ends at 31 s, at the corner of the search box (y on its lower bound, depth on its
# upper edge), because the P energy of this earthquake only builds up about 30 s after the predicted
# arrival (the stack of its 58 P records stays below a fifth of its maximum until then); the best
# models of every N then kept a subevent at 28-31 s and the misfit rose from 1 to 2 subevents. In
# the other ten v3 runs the 1-subevent centroid lies in the first 8-92% of its range. "On the
# bound" is the last 1% of the range, as in plot_lcurve.py, which applies the same test to the best
# model of every N. The range is doubled, because the magnitude gave no better
# scale, and the P, SH and regional windows and the moment-rate window grow by the same number of
# seconds (the records are far longer than the windows; the body-wave windows stop at 300 s as in
# calibrate_from_1sub.sh, which then shortens range and windows again from the new 1-subevent result).
widen_centroid_range_while_on_bound() {
    local inversion_dir=$1 event=$2 attempt centroid upper lower added folder
    for attempt in 1 2 3; do
        centroid=$(gawk '{print $4}' "$inversion_dir/best_model_hybrid.dat")
        upper=$(gawk '$1=="cen_max"{print $2}' "$inversion_dir/search_par.file")
        lower=$(gawk '$1=="cen_min"{print $2}' "$inversion_dir/search_par.file")
        if [ "$(gawk -v c="$centroid" -v u="$upper" -v l="$lower" 'BEGIN{print (u - c < 0.01 * (u - l)) ? 1 : 0}')" -ne 1 ]; then
            return 0
        fi
        if [ "$attempt" -eq 3 ]; then
            echo "WARNING: the 1-subevent centroid ($centroid s) is still on the upper bound of its range ($upper s) after two doublings: check the origin time of the event, or set cen_max in the search_par.file of every inv_${event}_*sub folder by hand"
            return 0
        fi
        added=$(gawk -v u="$upper" -v l="$lower" 'BEGIN{printf "%d", u - l + 0.999}')
        echo "centroid range: the 1-subevent centroid ($centroid s) is on the upper bound of its range ($upper s); range and windows of every N lengthened by $added s, and N = 1 runs again"
        for folder in inv_"${event}"_*sub; do
            [ -s "$folder/search_par.file" ] && [ -s "$folder/Par.file" ] || continue
            gawk -v a="$added" '$1=="cen_max"{$2 = $2 + a} {print}' "$folder/search_par.file" > "$folder/search_par.file.tmp" \
                && mv "$folder/search_par.file.tmp" "$folder/search_par.file" || return 1
            # the regional window keeps the cap of gen_par.py (500 samples from its start)
            gawk -v a="$added" '
                NR == FNR { if ($1 == "bg_timeRayl=") regional_begin = $2; if ($1 == "deltat=") interval = $2; next }
                /^(nd_timeP|nd_timeSH|nd_timeRayl|stf_etime)=/ {
                    key = $1; value = $2 + a
                    if ((key == "nd_timeP=" || key == "nd_timeSH=") && value > 300) value = 300
                    if (key == "nd_timeRayl=" && value > regional_begin + int(499 * interval)) value = regional_begin + int(499 * interval)
                    comment = $0; sub(/^[^#]*/, "", comment); print key " " value " " comment; next }
                { print }' "$folder/Par.file" "$folder/Par.file" > "$folder/Par.file.tmp" && mv "$folder/Par.file.tmp" "$folder/Par.file" || return 1
        done
        run_inversion "$inversion_dir" || return 1
    done
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
        # warm start: half of the chains begin from the best N-1 model plus one subevent (finv.f90);
        # each N refines the previous one instead of searching from scratch
        rm -f "$inversion_dir/warm_start.dat"
        previous_best="inv_${event_name}_$((subevent_count - 1))sub/best_model_hybrid.dat"
        [ "$subevent_count" -gt 1 ] && [ -s "$previous_best" ] && cp "$previous_best" "$inversion_dir/warm_start.dat"
        # a failed N stops step 3: the L-curve would otherwise select among the N that happened to work
        run_inversion "$inversion_dir" || { echo "step3 ABORTED at $inversion_dir"; exit 1; }
        if [ -z "$requested_counts" ] && [ "$subevent_count" -eq 1 ] && [ -f "inv_${event_name}_2sub/Par.file" ]; then
            widen_centroid_range_while_on_bound "$inversion_dir" "$event_name" \
                || { echo "step3 ABORTED: widening the centroid range failed for $event_name"; exit 1; }
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
