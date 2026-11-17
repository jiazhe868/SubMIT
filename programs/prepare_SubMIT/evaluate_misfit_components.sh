#!/bin/bash
# evaluate_misfit_components.sh INVERSION_DIR - the penalized misfit of the best sample of every chain,
# split into the data misfit and the two penalty factors: INVERSION_DIR/misfit_components.txt, which
# plot_lcurve.py tabulates in lcurve.txt and draws next to the L-curve (README section 4.6).
#
# Why. finv samples one number per model,
#     penalized misfit = data misfit x exp(weighted CLVD / DCconstrain) x exp(stress penalty)
# (forward_model.c, README section 4.4), so the L-curve alone cannot say whether N + 1 subevents fit
# the data better or only pay a smaller penalty. Split this way, the v3 reruns of the four examples
# (2026-11-16) showed three things. The stress factor was 1.000 for every best model, while the worst
# P or T axis sat at 26-30 deg from 2 subevents on in California and Venezuela and from 4 on in Chile:
# a limit the models lean on, not a cost they pay. California's choice
# between 4 and 5 subevents rested on the data misfit alone (two searches 5.4% apart at N = 5, their
# CLVD factors within 0.2%). Chile's two 3-subevent searches ended in different trade-offs (data
# misfit 0.697 x CLVD factor 1.150 against 0.835 x 1.005).
#
# How. ffwd re-evaluates the first line of every NNNbest.dat (the chain's best sample; finv writes it
# best-first) after one read of the Green's functions (about 20 s; each model then takes well under
# a second). A folder without chains (figures mode of an example ships them for the selected N only)
# has best_model_hybrid.dat alone. A run that was sampled with other code re-evaluates differently
# for every chain (California's v2 samples under the v3 timing: 5%); its factors would not multiply
# to the L-curve, so no table is written when the median chain is more than 1% off. The test is on
# the median because a single stored sample can re-evaluate a little off: models are written with 4
# decimals, and next to a jump of the misfit (an alignment lag that changes by one sample) that
# rounding moved one of the 32 chain bests of the 32-rank Venezuela 1-subevent run by 0.145%, while
# the other 31 agreed within 0.001%. The table lists both misfits of every chain, and the summary
# line counts the chains that are more than 0.1% off. Failures only warn: the table is a diagnostic,
# and step 3 goes on.
#
# usage (in IRIS/; step3_do_inversions.sh calls it after each inversion, with the compiler
# environment of programs/compiler_environment.sh loaded):
#     bash ../programs/prepare_SubMIT/evaluate_misfit_components.sh inv_<event>_<N>sub
# FFWD_PROGRAM=path/to/ffwd uses that binary instead of building ffwd in the inversion folder, for a
# run whose ffwd.f90 predates the model lists (2026-11-16).
set -u
inversion_dir=${1:?usage: evaluate_misfit_components.sh INVERSION_DIR}
ffwd_program=""
if [ -n "${FFWD_PROGRAM:-}" ]; then ffwd_program=$(readlink -f "$FFWD_PROGRAM"); fi
cd "$inversion_dir" || { echo "WARNING: misfit components skipped: no folder $inversion_dir"; exit 0; }
components_file=misfit_components.txt
models_file=misfit_components_models.txt
sampled_file=misfit_components_sampled.txt
log_file=misfit_components.log
rm -f "$components_file" "$models_file" "$sampled_file" "$log_file"

# ---------------------------------------------------------------------------------------------------
# the models: the best sample of every chain, in chain order (Input.model format, blank-line separated)
# ---------------------------------------------------------------------------------------------------
best_files=$(ls [0-9]*best.dat 2> /dev/null | sort -n)
if [ -z "$best_files" ] && [ -s best_model_hybrid.dat ]; then best_files=best_model_hybrid.dat; fi
[ -n "$best_files" ] || { echo "WARNING: misfit components skipped: no best sample in $inversion_dir"; exit 0; }
# a best-sample line: subevent count, "<neq", then 8 numbers per subevent, the first the misfit and
# the other 7 the Input.model columns (the same mapping as step4_figures.sh)
gawk -v models_file="$models_file" -v sampled_file="$sampled_file" 'FNR == 1 {
    chain = (FILENAME == "best_model_hybrid.dat") ? "best" : FILENAME
    sub(/best\.dat$/, "", chain)
    for (subevent = 0; subevent < $1; subevent++) {
        first = 3 + 8 * subevent
        print $(first + 1), $(first + 2), $(first + 3), $(first + 4), $(first + 5), $(first + 6), $(first + 7) > models_file
    }
    print "" > models_file
    print chain, $3 > sampled_file
    nextfile
}' $best_files

# ---------------------------------------------------------------------------------------------------
# one ffwd run for all models (compact: only the misfit lines, no files written)
# ---------------------------------------------------------------------------------------------------
if [ -z "$ffwd_program" ]; then
    # the objects of finv's build are reused: only ffwd.f90 is compiled
    make ffwd > make_ffwd.log 2>&1 || { echo "WARNING: misfit components skipped: ffwd did not build in $inversion_dir (see make_ffwd.log)"; exit 0; }
    ffwd_program=./ffwd
fi
FFWD_COMPACT=1 "$ffwd_program" "$models_file" > "$log_file" 2>&1 < /dev/null \
    || { echo "WARNING: misfit components skipped: ffwd failed in $inversion_dir (see $log_file)"; exit 0; }

# ---------------------------------------------------------------------------------------------------
# the table: per model the "bandresid:" line (residual and norm sums of P, SH, regional) and the
# "misfit_components:" line (data misfit, CLVD factor, stress factor, product, worst P and T axes)
# ---------------------------------------------------------------------------------------------------
gawk -v sampled_file="$sampled_file" -v components_file="$components_file" -v folder="$inversion_dir" '
    function ratio_or_nan(residual, norm) { return (norm > 0) ? sprintf("%.4f", residual / norm) : "nan" }
    BEGIN {
        while ((getline sampled_line < sampled_file) > 0) {
            split(sampled_line, sampled_fields, " ")
            model_count++
            chain[model_count] = sampled_fields[1]
            sampled_misfit[model_count] = sampled_fields[2]
        }
    }
    /^bandresid:/ { for (wave = 0; wave < 3; wave++) { residual[wave] = $(2 + wave); norm[wave] = $(6 + wave) } next }
    /^misfit_components:/ {
        evaluated++
        table_line[evaluated] = sprintf("%s %s %.6f %.6f %.6f %.6f %.2f %.2f %s %s %s", chain[evaluated],
            sampled_misfit[evaluated], $5, $2, $3, $4, $6, $7, ratio_or_nan(residual[0], norm[0]),
            ratio_or_nan(residual[1], norm[1]), ratio_or_nan(residual[2], norm[2]))
        deviation[evaluated] = $5 / sampled_misfit[evaluated] - 1
        if (deviation[evaluated] < 0) deviation[evaluated] = -deviation[evaluated]
        if (deviation[evaluated] > 1e-3) chains_off_count++
        if (best == 0 || sampled_misfit[evaluated] + 0 < sampled_misfit[best] + 0) {
            best = evaluated; best_summary = sprintf("data %.4f x CLVD %.4f x stress %.4f = %.4f", $2, $3, $4, $5)
        }
    }
    END {
        if (evaluated != model_count) {
            printf "WARNING: misfit components skipped: ffwd evaluated %d of %d models in %s\n", evaluated, model_count, folder
            exit 1
        }
        asort(deviation, sorted_deviation)
        median_deviation = sorted_deviation[int((evaluated + 1) / 2)]
        if (median_deviation > 1e-2) {
            printf "WARNING: misfit components skipped: the re-evaluated misfits in %s differ from the sampled ones by %.2f%% for the median chain (sampled with other code?)\n", folder, 100 * median_deviation
            exit 1
        }
        print "chain sampled_misfit penalized_misfit data_misfit clvd_factor stress_factor worst_p_deg worst_t_deg p_misfit sh_misfit regional_misfit" > components_file
        for (index_ = 1; index_ <= evaluated; index_++) print table_line[index_] > components_file
        printf "misfit components: %s, %d chains; best %s", folder, evaluated, best_summary
        if (chains_off_count > 0) printf "; %d chain(s) re-evaluated more than 0.1%% from the sampled misfit", chains_off_count
        printf "\n"
    }' "$log_file"
rm -f "$models_file" "$sampled_file"
exit 0
