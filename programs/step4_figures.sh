#!/bin/bash
# step4_figures.sh N - forward model, posterior summary and figures of the N-subevent result.
#
# For every event in event_list.dat: the best sample of inv_<event>_<N>sub becomes Input.model;
# hist/posterior_summary.py pools the converged chains (allsamples.dat, histograms, misfit
# evolution, convergence.txt); fwd_<event>_<N>sub is rebuilt from the inversion folder, so the
# forward model sees exactly the data and configuration that were sampled; ffwd writes waveforms,
# moment tensors and per-station residuals; the matplotlib scripts draw the figures; and
# figs_and_results/ collects them with RESULTS.md (set SUBMIT_NO_PUBLISH=1 to regenerate a
# non-selected N without overwriting figs_and_results/).
#
# History: replaces step4_fwd_tempnsub.sh (v2), which overlaid the step-2 forward template with two
# source trees (a stale parser once segfaulted on newer Par.file keys) and kept a GMT4 plotting path.
#
# usage (in IRIS/):  bash ../programs/step4_figures.sh 4
if [ -z "$BASH_VERSION" ]; then exec bash "$0" "$@"; fi
set -u
subevent_count=${1:?usage: step4_figures.sh N}
set --
. ../programs/compiler_environment.sh
python_executable=$SUBMIT_PYTHON
case "$python_executable" in */intel*|*oneapi*) python_executable=/usr/bin/python3 ;; esac
products_missing=0

while read -r event_name rest; do
    inversion_dir="inv_${event_name}_${subevent_count}sub"
    forward_dir="fwd_${event_name}_${subevent_count}sub"
    best_model="$inversion_dir/best_model_hybrid.dat"
    if [ ! -s "$best_model" ]; then
        sort -gk3 "$inversion_dir"/[0-9]*best.dat 2>/dev/null | head -1 > "$best_model"
    fi
    [ -s "$best_model" ] || { echo "step4: no inversion result in $inversion_dir"; products_missing=1; continue; }
    # best sample -> Input.model (one row per subevent: cen x y dur vr theta depth) + misfit.dat
    gawk '{for (i = 0; i < $1; i++) {b = 3 + 8*i; print $(b+1), $(b+2), $(b+3), $(b+4), $(b+5), $(b+6), $(b+7)}}' "$best_model" > "$inversion_dir/Input.model"
    gawk '{print $3}' "$best_model" > "$inversion_dir/misfit.dat"

    ( cd "$inversion_dir" && "$python_executable" ../../programs/hist/posterior_summary.py > posterior_summary.log 2>&1 ) \
        || { echo "step4: posterior summary failed (see $inversion_dir/posterior_summary.log)"; products_missing=1; }

    rm -rf "$forward_dir"
    mkdir "$forward_dir"
    ( cd "$inversion_dir" && cp -r $(ls | grep -v -e 'chain\.dat$' -e '^pooled_ensemble\.dat$' -e '\.o$' -e '^finv$' -e '^mpi\.out$') "../$forward_dir/" )
    cp ../programs/code_SubMIT/{forward_model.c,ffwd.f90,finv.f90,mtdcmp.f,sacsubc.c,sacsubc.h,Makefile} "$forward_dir/"
    ( cd "$forward_dir" && make clean > /dev/null 2>&1 && make ffwd > make.log 2>&1 ) \
        || { echo "step4: COMPILE FAILED in $forward_dir (see make.log)"; products_missing=1; continue; }
    ( cd "$forward_dir" && ./ffwd > ffwd.log 2>&1 < /dev/null \
        && "$python_executable" ../../programs/fwd_SubMIT/plot_waveform_fits_mpl.py > figures.log 2>&1 \
        && "$python_executable" ../../programs/fwd_SubMIT/plot_subevent_map_mpl.py >> figures.log 2>&1 \
        && "$python_executable" ../../programs/fwd_SubMIT/plot_station_map_mpl.py >> figures.log 2>&1 ) \
        || echo "step4: forward run or figures failed in $forward_dir (see ffwd.log, figures.log)"
    for product in "$forward_dir/fits_P.pdf" "$forward_dir/fits_SH.pdf" "$forward_dir/subevents.pdf" \
                   "$inversion_dir/histoplot_py.pdf" "$forward_dir/fm.dat"; do
        [ -s "$product" ] || { echo "MISSING step4 product: $product"; products_missing=1; }
    done

    if [ "${SUBMIT_NO_PUBLISH:-0}" != "1" ]; then
        mkdir -p figs_and_results
        cp "$forward_dir"/fits_*.pdf "$forward_dir"/subevents.pdf "$forward_dir"/subevents.png \
           "$forward_dir"/stations.png "$forward_dir"/stations.pdf figs_and_results/ 2>/dev/null
        cp "$inversion_dir"/histoplot_py.pdf "$inversion_dir"/misfit_evolution.pdf "$inversion_dir"/convergence.txt figs_and_results/ 2>/dev/null
        cp lcurve.pdf lcurve.txt figs_and_results/ 2>/dev/null
        "$python_executable" ../programs/fwd_SubMIT/make_results_doc.py "$event_name" "$subevent_count" \
            && echo "results document -> figs_and_results/RESULTS.md"
    fi
done < event_list.dat
[ "$products_missing" -eq 0 ] || { echo "step4 FAILED: products missing (see above)"; exit 1; }
echo "step4 done"
