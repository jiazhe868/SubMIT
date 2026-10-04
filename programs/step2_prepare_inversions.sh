#!/bin/bash
# step2: prepare inversions. Fails fast with named errors instead of cascading.
# The manual invokes this as "sh step2...", which may be dash: re-exec under bash.
if [ -z "$BASH_VERSION" ]; then exec bash "$0" "$@"; fi
set -o pipefail
exec > >(tee step2.log) 2>&1

echo "== step2 preflight =="
fail=0
for tool in sac saclst gawk xargs; do
    command -v "$tool" > /dev/null || { echo "MISSING tool: $tool"; fail=1; }
done
for bin in ../programs/fk3.2/fk ../programs/mtel3/mtel3; do
    [ -x "$bin" ] || { echo "MISSING/uncompiled binary: $bin (run ./submit build in the repository root)"; fail=1; }
done
python -c "import numpy, scipy, pandas, requests, obspy" 2>/dev/null \
    || { echo "MISSING python packages (need numpy scipy pandas requests obspy; run step2 in the conda env)"; fail=1; }
ls -1 -d [12]???-??-??* > /dev/null 2>&1 || { echo "NO event folders ([12]???-??-??*) here"; fail=1; }
for d in $(ls -1 -d [12]???-??-??* 2>/dev/null); do
    [ -d "../IRISloc/${d}_loc/data" ] || echo "WARNING: no local data ../IRISloc/${d}_loc/data for $d"
done
[ "$fail" -eq 0 ] || { echo "== preflight FAILED, aborting =="; exit 1; }
echo "== preflight OK =="

stage() { echo; echo "== step2 stage: $* =="; }

stage "event list"
ls -1 -d [12]???-??-??* | gawk '{print "sh ../programs/prepare_SubMIT/gen_event_list.sh "$1 }'  | sh > event_list.dat
[ -s event_list.dat ] || { echo "ERROR: event_list.dat empty"; exit 1; }
stage "copy code + data"
cat event_list.dat  | gawk '{print "rm -rf inv_"$1"* fwd_"$1"*";print "cp -r ../programs/code_SubMIT ./inv_"$1; print "cp -r "$1"/data/ inv_"$1"/"; print "cp -r ../IRISloc/"$1"_loc/data/ inv_"$1"/dataloc"}' | sh
stage "integrate teleseismic data"
cat event_list.dat | gawk '{print "sh ../programs/prepare_SubMIT/process_tel_comp.sh inv_"$1}' | sh
stage "Green's functions (this is the long stage)"
cat event_list.dat  | gawk '{print "rm -rf gf_"$1; print "cp -r ../programs/gf_SubMIT gf_"$1; print "cd gf_"$1"/loc";print "sh doloc.sh";print "cd ../";print "sh dotel.sh";print "cd .."}' | sh
while read -r ev rest; do
    # find, not a bare glob: at wide depth apertures (>50 vmodel_* dirs) the
    # glob expansion exceeds ARG_MAX and ls fails on E2BIG despite complete GFs
    find gf_${ev}/greenFuncDir_disp -maxdepth 2 -name "????*.grn.????" -print -quit 2>/dev/null | grep -q . \
        || { echo "ERROR: no teleseismic Green's functions for $ev (gf_${ev}/greenFuncDir_disp)"; exit 1; }
    # local GFs: every search depth must have transferred double-couple files
    for d in gf_${ev}/greenFuncDir_disp/vmodel_*; do
        ls "$d"/[0-9]*.grn.ddpz > /dev/null 2>&1 \
            || { echo "ERROR: no local Green's functions in $d (fk failed? see log)"; exit 1; }
    done
done < event_list.dat
stage "station screening + station info"
cat event_list.dat | gawk '{print "cd inv_"$1;print "python ../../programs/prepare_SubMIT/screen_station_amplitudes.py";print "cp ../gf_"$1"/loc/distdep.dat .";print "sh prepare_stationinfo.sh";print "cd .."}' | sh
stage "aftershock prior + weights + Par files"
cat event_list.dat | gawk '{print "cp ../programs/prepare_SubMIT/*.py inv_"$1; print "cd inv_"$1;print "python fetch_aftershocks.py "$2"T"$3,$5,$6,$7,$4; print "python caldens.py "$5,$6; print "python gen_weight.py "$4; print "python gen_par.py "$1,$4,$5,$6,$7; print "python gen_totalmt.py "$2"T"$3,$5,$6,$7,$4; print "cd ../";}' | sh
while read -r ev rest; do
    for f in edges.dat seisdens.dat weights.dat Par.file; do
        [ -s "inv_${ev}/$f" ] || { echo "ERROR: inv_${ev}/$f missing/empty (prior or Par stage failed)"; exit 1; }
    done
done < event_list.dat
stage "replicate 1-5 subevent folders + search_par"
while read line; do
    event_name=$(echo $line | awk '{print $1}')
    mw=$(echo $line | awk '{print $4}')
    nsub=2
    while [ $nsub -le 5 ]; do
        cp -r "inv_${event_name}" "inv_${event_name}_${nsub}sub"
        nsub=$((nsub + 1))
    done
    cp -r "inv_${event_name}" "fwd_${event_name}"
    mv "inv_${event_name}" "inv_${event_name}_1sub"
    nsub=1
    while [ $nsub -le 5 ]; do
        cd "inv_${event_name}_${nsub}sub"
        python get_search_par.py "$event_name" "$mw" "$nsub"
        cd ..
        nsub=$((nsub + 1))
    done
done < event_list.dat

#source /opt/intel/oneapi/setvars.sh
