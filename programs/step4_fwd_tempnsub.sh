#!/bin/bash

if [ -z "$1" ]; then
    echo "No argument provided. Exiting."
    exit 1
else
    nn=$1
fi
nn=$1

cat event_list.dat | gawk -v nn="$nn" '{print "cp ../programs/fwd_SubMIT/get_best_model.py inv_"$1"_"nn"sub"; print "cd inv_"$1"_"nn"sub"; print "python get_best_model.py"; print "cd .."}' | sh

# remember the user's (conda) python BEFORE the Intel env hijacks PATH with
# intelpython (broken numpy/matplotlib ABI) - the histogram step needs it
PYEXE=${SUBMIT_PYTHON:-$(command -v python3 || command -v python)}
case "$PYEXE" in */intel*|*oneapi*) PYEXE=/usr/bin/python3 ;; esac

# validated Intel env (same logic as step3); generic MPI wrappers work too
# SUBMIT_NO_INTEL=1 skips this (e.g. when building with GNU mpicc/mpif90:
# an Intel mpirun cannot launch an OpenMPI/MPICH-built binary)
[ "${SUBMIT_NO_INTEL:-0}" = "1" ] && INTEL_CANDIDATES="" || \
    INTEL_CANDIDATES="$INTEL_SETVARS /opt/intel/oneapi/setvars.sh $HOME/intel/oneapi/setvars.sh"
for sv in $INTEL_CANDIDATES; do
    [ -n "$sv" ] && [ -f "$sv" ] || continue
    . "$sv" --force > /dev/null 2>&1
    command -v ifx > /dev/null 2>&1 || for cv in $(ls -1d "$(dirname "$sv")"/compiler/*/env/vars.sh 2>/dev/null | sort -rV); do
        [ -f "$cv" ] && . "$cv" > /dev/null 2>&1 && command -v ifx > /dev/null 2>&1 && break
    done
    command -v ifx > /dev/null 2>&1 && break
done

# plotting: default is the GMT-free matplotlib path (portable - python+obspy
# are already required by steps 1-2). Set SUBMIT_PLOTTER=gmt for the legacy
# GMT4+pssac2 scripts, which then get a named preflight.
SUBMIT_PLOTTER=${SUBMIT_PLOTTER:-mpl}
if [ "$SUBMIT_PLOTTER" = "gmt" ]; then
    pfail=0
    for t in pssac2 psbasemap pstext psmeca psxy minmax gmtset; do
        command -v "$t" > /dev/null || { echo "MISSING plot tool: $t (GMT4 + pssac2 required for SUBMIT_PLOTTER=gmt)"; pfail=1; }
    done
    [ "$pfail" -eq 0 ] || { echo "step4: GMT plotting toolchain incomplete - aborting"; exit 1; }
fi

pids=""
# Read each line from event_list.dat, extract the first field ($1), and process in parallel
# no pipeline here: a pipeline runs the loop in a SUBSHELL, orphaning the
# backgrounded blocks so that "wait" (and the manifest) cannot see them
while read -r event_id rest; do
  {
    # Per-nsub fwd directory: fwd_<event> stays a pristine step2 template and
    # different subevent counts no longer overwrite each other's figures
    fdir="fwd_${event_id}_${nn}sub"
    rm -rf "$fdir"
    cp -r "fwd_$event_id" "$fdir"
    cp -r ../programs/fwd_SubMIT/* "$fdir"
    # always ship the CURRENT shared sources (a stale positional-parser
    # sub_init.c in the step2-era template segfaulted on newer Par.file keys)
    cp ../programs/code_SubMIT/sub_init.c ../programs/code_SubMIT/sub_header.h "$fdir"
    cp "inv_${event_id}_${nn}sub/Input.model" "inv_${event_id}_${nn}sub/misfit.dat" "inv_${event_id}_${nn}sub/search_par.file" "inv_${event_id}_${nn}sub/stations.info" "inv_${event_id}_${nn}sub/stationsSH.info" "inv_${event_id}_${nn}sub/stationsloc.info"  "inv_${event_id}_${nn}sub/Par.file" "$fdir"
    cp "inv_${event_id}_${nn}sub/mainshock.dat" "$fdir" 2>/dev/null
    # data MUST come from the inv dir, not the step2-era template: a stale
    # dataloc (e.g. 400-s vs 600-s traces) makes ffwd re-solve DIFFERENT
    # moment tensors than the inversion sampled - the plotted mechanisms
    # would not be the inverted ones
    rm -rf "$fdir/data" "$fdir/dataloc"
    cp -r "inv_${event_id}_${nn}sub/data" "inv_${event_id}_${nn}sub/dataloc" "$fdir/"
    cp "inv_${event_id}_${nn}sub/totalmt.dat" "$fdir/" 2>/dev/null
    # regenerate the convergence-gated pool from THIS run's chains (a stale
    # pooled_ensemble.dat from an earlier run skews the error bars), then
    # histograms: allsamples.dat feeds the 95% error bars on the maps
    cp ../programs/prepare_SubMIT/pool_chains.py ../programs/hist/process_and_plot.py \
        ../programs/hist/plot_misfit_evolution.py \
        "inv_${event_id}_${nn}sub/" || { echo "MISSING hist-stage scripts"; exit 1; }
    (cd "inv_${event_id}_${nn}sub" && "$PYEXE" pool_chains.py > pool.log 2>&1 \
        && "$PYEXE" process_and_plot.py > hist.log 2>&1 \
        && "$PYEXE" plot_misfit_evolution.py >> hist.log 2>&1) \
        || { echo "hist stage FAILED (see inv_${event_id}_${nn}sub/pool.log, hist.log)"; exit 1; }
    cp "inv_${event_id}_${nn}sub/allsamples.dat" "$fdir" 2>/dev/null
    # Navigate to the fwd directory and run commands
    cd "$fdir" || exit
    rm -f waveforms/*.sac 2>/dev/null   # no stale traces from previous runs
    # make clean first: copied stale .o files can tie-timestamp with fresh
    # sources and silently link an outdated parser into ffwd
    make clean > /dev/null 2>&1
    make || { echo "COMPILE FAILED in $fdir"; exit 1; }
    [ -x ./ffwd ] || { echo "ffwd missing in fwd_$event_id"; exit 1; }
    ./ffwd
    
    # Figures: matplotlib (default, portable) or legacy GMT4 (SUBMIT_PLOTTER=gmt)
    if [ "$SUBMIT_PLOTTER" = "gmt" ]; then
        cd plotP/ && sh plotP.sh
        cd ../plotPvel/ && sh plotPvel.sh
        cd ../plotSH/ && sh plotSH.sh
        cd ../plotrayl/ && sh plotrayl.sh
        cd ../subeveplot/ && sh plot.sh
        cd ../../
    else
        "$PYEXE" plot_waveform_fits_mpl.py || echo "WARNING: waveform-fit plates failed"
        "$PYEXE" plot_subevent_map_mpl.py || echo "WARNING: subevent figure failed"
        "$PYEXE" plot_station_map_mpl.py || echo "WARNING: station map failed"
        cd ../
    fi

  } > output_${event_id}.dat &  # Run in background (child of THIS shell)
done < event_list.dat

# Wait for all background processes to complete. NOTE: the per-event blocks are
# children of the "awk | while" pipeline SUBSHELL, so $pids is empty out here
# (the old loop waited on nothing and the driver exited while plots were still
# rendering). A bare "wait" waits on ALL child processes of this shell.
wait

# Product manifest: the per-event work runs in background subshells whose exit
# codes are unobservable here, so verify the artifacts themselves.
mfail=0
while read -r event_id rest; do
    fdir="fwd_${event_id}_${nn}sub"
    if [ "$SUBMIT_PLOTTER" = "gmt" ]; then
        want="$fdir/plotP/map.ps $fdir/plotSH/map.ps $fdir/subeveplot/fig.ps"
    else
        want="$fdir/fits_P.pdf $fdir/fits_SH.pdf $fdir/fits_rayl.pdf $fdir/subevents.pdf"
    fi
    for f in $want "inv_${event_id}_${nn}sub/histoplot_py.pdf"; do
        [ -s "$f" ] || { echo "MISSING step4 product: $f"; mfail=1; }
    done
    nw=$(ls "$fdir"/waveforms/*.sac 2>/dev/null | wc -l)
    [ "$nw" -gt 0 ] || { echo "MISSING step4 product: $fdir/waveforms/*.sac"; mfail=1; }
done < event_list.dat
[ "$mfail" -eq 0 ] || { echo "step4 FAILED: products missing (see above)"; exit 1; }
echo "step4 product manifest: OK"

# collect the selected-nsub figures + a formatted results document
# SUBMIT_NO_PUBLISH=1 skips this block: publishing is per-nn, so regenerating
# a NON-selected nn (e.g. refreshing stale fwd dirs) would otherwise clobber
# figs_and_results with the wrong subevent count
if [ "${SUBMIT_NO_PUBLISH:-0}" = "1" ]; then
    echo "SUBMIT_NO_PUBLISH=1: skipping figs_and_results publish"
    echo "All processes completed."
    exit 0
fi
while read -r event_id rest; do
    fdir="fwd_${event_id}_${nn}sub"
    out="figs_and_results"
    mkdir -p "$out"
    cp "$fdir"/fits_*.pdf "$fdir"/subevents.pdf "$fdir"/subevents.png "$out"/ 2>/dev/null
    cp "$fdir"/stations.png "$fdir"/stations.pdf "$out"/ 2>/dev/null
    cp "inv_${event_id}_${nn}sub/histoplot_py.pdf" \
       "inv_${event_id}_${nn}sub/misfit_evolution.pdf" "$out"/ 2>/dev/null
    cp lcurve.pdf lcurve.txt "$out"/ 2>/dev/null
    "$PYEXE" ../programs/fwd_SubMIT/make_results_doc.py "$event_id" "$nn" && \
        echo "results document -> $out/RESULTS.md"
done < event_list.dat
echo "All processes completed."

