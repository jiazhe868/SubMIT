#!/bin/bash
# step2_prepare_inversions.sh - everything between the processed records and the MCMC runs.
#
# For every event folder in IRIS/ (output of step 1): event list, a copy of the inversion code and
# the data per event, integration of the teleseismic records to displacement, teleseismic and
# regional Green's functions, station screening (amplitude, signal-to-noise, drift, clipping, dense
# networks thinned to 100 stations) and station lists, aftershock prior, wave-type weights, Par.file,
# total moment tensor, and one inversion folder per number of subevents (1-5) with its
# search_par.file. Each stage fails loudly with a named error instead of cascading.
#
# usage (in IRIS/):  bash ../programs/step2_prepare_inversions.sh
if [ -z "$BASH_VERSION" ]; then exec bash "$0" "$@"; fi
set -o pipefail
exec > >(tee step2.log) 2>&1

echo "== step2 preflight =="
preflight_failed=0
for tool in sac saclst gawk xargs; do
    command -v "$tool" > /dev/null || { echo "MISSING tool: $tool"; preflight_failed=1; }
done
for program in ../programs/fk3.2/fk ../programs/mtel3/mtel3; do
    [ -x "$program" ] || { echo "MISSING/uncompiled binary: $program (run ./submit build in the repository root)"; preflight_failed=1; }
done
python -c "import numpy, scipy, pandas, requests, obspy" 2>/dev/null \
    || { echo "MISSING python packages (need numpy scipy pandas requests obspy; run step2 in the conda env)"; preflight_failed=1; }
ls -1 -d [12]???-??-??* > /dev/null 2>&1 || { echo "NO event folders ([12]???-??-??*) here"; preflight_failed=1; }
for event_folder in $(ls -1 -d [12]???-??-??* 2>/dev/null); do
    [ -d "../IRISloc/${event_folder}_loc/data" ] || echo "WARNING: no regional data ../IRISloc/${event_folder}_loc/data for $event_folder"
done
[ "$preflight_failed" -eq 0 ] || { echo "== preflight FAILED, aborting =="; exit 1; }
echo "== preflight OK =="

stage() { echo; echo "== step2 stage: $* =="; }

stage "event list"
for event_folder in $(ls -1 -d [12]???-??-??*); do
    sh ../programs/prepare_SubMIT/gen_event_list.sh "$event_folder"
done > event_list.dat
[ -s event_list.dat ] || { echo "ERROR: event_list.dat empty"; exit 1; }

while read -r event_name origin_date origin_time magnitude longitude latitude depth; do
    inversion_dir="inv_${event_name}"
    stage "$event_name: inversion code + data"
    rm -rf "inv_${event_name}"* "fwd_${event_name}"*
    cp -r ../programs/code_SubMIT "$inversion_dir"
    cp -r "${event_name}/data/" "$inversion_dir/"
    cp -r "../IRISloc/${event_name}_loc/data/" "$inversion_dir/dataloc"

    stage "$event_name: teleseismic records integrated to displacement"
    sh ../programs/prepare_SubMIT/process_tel_comp.sh "$inversion_dir" > "$inversion_dir/process_tel_comp.log" 2>&1 \
        || { echo "ERROR: integration to displacement failed for $event_name (see $inversion_dir/process_tel_comp.log)"; exit 1; }

    stage "$event_name: Green's functions (the long stage)"
    rm -rf "gf_${event_name}"
    cp -r ../programs/gf_SubMIT "gf_${event_name}"
    # bash, not sh: both scripts cap their parallel jobs with "jobs -rp", which dash (sh on Debian and
    # Ubuntu) rejects ("Illegal option -r", 27 to 168 times in every run log until 2026-11-17), so the
    # cap never held and every depth of the grid was started at once
    ( cd "gf_${event_name}/loc" && bash doloc.sh ) || { echo "ERROR: regional Green's functions failed for $event_name"; exit 1; }
    ( cd "gf_${event_name}" && bash dotel.sh ) || { echo "ERROR: teleseismic Green's functions failed for $event_name"; exit 1; }
    # find, not a glob: at wide depth apertures (> 50 depth folders) the glob exceeds ARG_MAX
    find "gf_${event_name}/greenFuncDir_disp" -maxdepth 2 -name "????*.grn.????" -print -quit 2>/dev/null | grep -q . \
        || { echo "ERROR: no teleseismic Green's functions for $event_name (gf_${event_name}/greenFuncDir_disp)"; exit 1; }
    for depth_folder in "gf_${event_name}"/greenFuncDir_disp/vmodel_*; do
        ls "$depth_folder"/[0-9]*.grn.ddpz > /dev/null 2>&1 \
            || { echo "ERROR: no regional Green's functions in $depth_folder (fk failed? see the log)"; exit 1; }
    done

    stage "$event_name: station screening + station lists"
    ( cd "$inversion_dir" && python ../../programs/prepare_SubMIT/screen_station_amplitudes.py \
        && cp "../gf_${event_name}/loc/distdep.dat" . && python prepare_stationinfo.py ) \
        || { echo "ERROR: station screening failed for $event_name"; exit 1; }

    stage "$event_name: aftershock prior, weights, Par.file, total moment tensor"
    cp ../programs/prepare_SubMIT/*.py "$inversion_dir/"
    iso_time="${origin_date}T${origin_time}"
    ( cd "$inversion_dir" \
        && python fetch_aftershocks.py "$iso_time" "$longitude" "$latitude" "$depth" "$magnitude" \
        && python caldens.py "$longitude" "$latitude" \
        && python gen_weight.py "$magnitude" \
        && python gen_par.py "$event_name" "$magnitude" "$longitude" "$latitude" "$depth" \
        && python gen_totalmt.py "$iso_time" "$longitude" "$latitude" "$depth" "$magnitude" ) \
        || { echo "ERROR: prior/weights/Par stage failed for $event_name"; exit 1; }
    for required in edges.dat seisdens.dat weights.dat Par.file totalmt.dat; do
        [ -s "$inversion_dir/$required" ] || { echo "ERROR: $inversion_dir/$required missing or empty"; exit 1; }
    done

    stage "$event_name: one inversion folder per number of subevents (1-5)"
    for subevent_count in 2 3 4 5; do
        cp -r "$inversion_dir" "${inversion_dir}_${subevent_count}sub"
    done
    cp -r "$inversion_dir" "fwd_${event_name}"          # forward template for the step-3 calibration
    mv "$inversion_dir" "${inversion_dir}_1sub"
    for subevent_count in 1 2 3 4 5; do
        ( cd "${inversion_dir}_${subevent_count}sub" && python get_search_par.py "$event_name" "$magnitude" "$subevent_count" ) \
            || { echo "ERROR: search_par.file not written for ${inversion_dir}_${subevent_count}sub"; exit 1; }
    done
done < event_list.dat
echo "step2 done"
