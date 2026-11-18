#!/usr/bin/env python3
"""L-curve (best penalized misfit against the number of subevents) and the selected number.

Rule (README section 4.6, since 2026-11-17): the Akaike information criterion of every N,
    AIC(N) = ln(misfit_N) / sigma^2 + 2 k N,
where sigma = 0.05 is the width of the posterior that finv.f90 samples (its sigma_relative) and
k = 10 counts the free parameters of a subevent: five moment-tensor components, centroid time,
duration, depth, x and y. (The first subevent of two or more has fewer: it sits at the hypocentre,
so only its centroid time, its duration and its tensor are free. That is the same for every N from
2 on, so one more subevent adds ten parameters; only the step from one subevent to two is charged
three more than it adds.) Every N whose AIC is within 2 of the lowest is tied with it, and the
selection is the smallest of them. In words: a subevent costs 2 k sigma^2 = 5% of misfit (exactly:
ln(misfit) has to fall by 0.05, a misfit 4.9% lower, to tie with one subevent fewer), and of the
numbers that the data cannot tell apart the simplest is taken.

Where the formula comes from. finv samples exp(-misfit / (2 sigma^2 misfit_min)): a model whose
misfit is 1% higher has e^-2 of the density, as for n = 1 / sigma^2 = 400 independent data whose
noise level is taken from the best fit. For least squares with an unknown noise level the criterion
is AIC = n ln(RSS) + 2 K (Akaike, 1974; Burnham & Anderson, 2002); with the misfit in the place of
the residual sum of squares and K = k N parameters this is the line above. The 5% is not a new
constant: it follows from the sigma of the posterior and the size of a subevent. The misfit is not
a plain sum of squares and sigma is an assumed noise level, so this is a consistent price for a
subevent, not an exact test.

Why not the earlier rule. Until 2026-11-16 the selection was the smallest N within 5% of the lowest
misfit of the whole grid. That rule charges nothing per subevent, and it refers every N to the
largest ones, whose searches are the hardest to finish. Measured on the runs of the examples:
- Venezuela selected 6 (published v2 run), 5 (v3, 16 ranks) and 6 (v3, 32 ranks). In the 32-rank run
  5 subevents were rejected because 7 fit 6.5% better, 3.3% per added subevent. The AIC selects 5 in
  all three runs.
- California selected 4 (v3, 16 ranks) and 5 (v3, 32 ranks, where 4 subevents missed the band by
  0.05%). The rule above selects 4 in both, and in a third search with other chain seeds.

Ties. Models whose AIC differs by less than 2 are equally supported (Burnham & Anderson, 2002).
lcurve.txt names every such N, and the smallest is selected: a larger N has to be clearly better,
not better by an amount that the next search takes back. Measured on California with 32 ranks and
two sets of chain seeds: the fifth subevent lowered the misfit by 4.8% in one search and by 5.0% in
the other, on the two sides of the 4.9% that it costs, so the lowest AIC was at 4 and then at 5.
Both searches name the tie of 4 and 5 and both select 4. The models say the same: the two
4-subevent models agree (subevents within 0.35 s and 9 km of each other), the two 5-subevent models
place their fourth subevent 46 km apart, and in a search with the centroid range widened from 16
to 22 s the fifth subevent raised the misfit by 1.3%. Venezuela's 6 subevents (0.2 above 5, 32 ranks) are the
other tie of the examples. The best misfit of an N is what a stochastic search found, an upper
bound of its minimum, and larger N are harder to search, so a tie can end in a deeper search.

Grid extension. ./submit adds one N while the lowest AIC sits at the largest N tried (not while the
selection does: a tie at the edge of the grid would otherwise stop the grid one step before an N
that is clearly better). The line "lowest AIC: N=" of lcurve.txt is what it reads.

A best model on the bound of its centroid range. An N whose best model has its last subevent in
the last 1% of the centroid range (cen_min to cen_max of search_par.file, the test with which step
3 widens the range of the 1-subevent run) is listed with "bound" in the place of its AIC and cannot
be selected; the lowest AIC and the ties are taken over the other N. Such a subevent is not a free
optimum: it fits later moment release, or arrivals that the model does not predict, from the edge
of what the range allows. Measured on the 2026 Vanuatu Mw 7.3 event, 121 km deep (v3, two sets of
chain seeds): three subevents fit 11% better than one, and in both searches the third sat at
28.0 s of a range that ends at 28 s, 53 km east in a box that ends at 53 km, at the shallow end of
the depth grid, where it lowered the misfit of the three regional stations and raised those of P
and SH. Both searches agreed on it, so the agreement of two searches does not catch this. Without
that N the selection is one subevent, and two subevents fit within 3.4% of it. (With the range
widened by hand to 55 s that subevent settles at 31.5 s on the top of the depth grid, and its P wave
arrives between the pP phases of the first two subevents: it looks like depth-phase energy that the
model does not predict. PITFALLS D5.) When every N is on
the bound nothing can be excluded, and a note says that the selection should not be trusted.
The 1% is where the best models of the v3 runs separate: those pressed against the end of their
range lie in its last 0.75% (Colombia before the range was widened, Vanuatu, and three California
models that are not selected), and the nearest other one lies 2.07% before the end (California,
15.70 s of a range from 1.5 to 16 s: a subevent that stays at 15.8 s when the range ends at 22 s).

Island bests. finv resamples its chains within islands of 8 and, in the last burn-in round, across
all chains; island_best_misfits.txt holds the best misfit of every island before that round. They
are listed, and drawn as open circles, to show how far independent searches of one N ended from each
other. They do not enter the rule. An island best can lie a little below the best misfit of its row
(second 32-rank California search, 6 subevents: 1.3504 against 1.3572): the best misfit is the best
sample after the burn-in, and the burn-in passes through models that the sampling does not visit
again. The v3 code of 2026-11-16 bootstrapped them (1000 resamples, a
winner below 80% flagged as unstable); that test was removed with the old rule. It flagged all three
32-rank examples (22.5%, 68% and 71%), because an island best is a snapshot at 75% of the burn-in and
not a finished search, and under the AIC it voted against the selection it was meant to test
(Venezuela: 6 subevents in 71% of the resamples, 5 selected from the finished chains).

Misfit components (2026-11-16). The misfit is a product,
    penalized misfit = data misfit x exp(weighted CLVD / DCconstrain) x exp(stress penalty)
(forward_model.c), so the L-curve alone cannot tell a better data fit from a smaller penalty. After
each inversion, step 3 re-evaluates the best sample of every chain and writes the three factors to
misfit_components.txt (evaluate_misfit_components.sh). This script lists them for the best model of
every N in lcurve.txt and draws them next to the L-curve, with the worst P/T-axis rotation that the
stress penalty acts on. In the v3 reruns of the four examples they showed that California's choice
between 4 and 5 subevents rested on the data misfit (two searches 5.4% apart at N = 5, their CLVD
factors within 0.2%), and that the stress factor was 1.000 for every best model while the worst
axis of most of them sat at 26-30 deg. Runs without the file (v2, or made before this date) get the
L-curve alone.

usage (in IRIS/, after step 3): python ../programs/prepare_SubMIT/plot_lcurve.py
writes lcurve.txt and lcurve.pdf; the last line of lcurve.txt names the selection.
"""
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

posterior_sigma = 0.05          # sigma_relative of finv.f90; tests/test_regressions.py checks that the two are equal
parameters_per_subevent = 10    # 5 moment-tensor components, centroid time, duration, depth, x, y
log_misfit_drop_per_subevent = 2 * parameters_per_subevent * posterior_sigma**2      # 0.05: what a subevent has to earn
tie_aic_difference = 2.0        # models closer than this are equally supported (Burnham & Anderson, 2002)
bound_share_of_range = 0.01     # "on the bound": the last 1% of the centroid range, as in step3_do_inversions.sh (a test checks that)
stress_cone_deg = 30.0          # forward_model.c: the stress penalty starts beyond 30 deg of axis rotation
factor_labels = {"data_misfit": "data misfit", "clvd_factor": "CLVD factor", "stress_factor": "stress factor",
                 "worst_axis_deg": "worst axis (°)"}


def chain_best_misfits(inversion_dir):
    """best misfit of every chain (first line of NNNbest.dat, finv writes it best-first). Folders
    that carry only best_model_hybrid.dat (the published examples ship chains for the selected N
    only) count as one chain."""
    values = []
    for best_file in glob.glob(os.path.join(inversion_dir, "[0-9]*best.dat")):
        fields = open(best_file).readline().split()
        if len(fields) > 2:
            values.append(float(fields[2]))
    snapshot = os.path.join(inversion_dir, "best_model_hybrid.dat")
    if not values and os.path.exists(snapshot):
        fields = open(snapshot).readline().split()
        if len(fields) > 2:
            values.append(float(fields[2]))
    return np.array(values)


def last_centroid_against_its_bound(inversion_dir):
    """{"last_centroid_s": centroid time of the last subevent of the folder's best model,
    "range_end_s": upper bound of the centroid range, "is_on_bound": whether that time lies in the
    last 1% of the range}. A folder without search_par.file or without a model (hand-made folders,
    the tests) is not on the bound."""
    unknown = {"last_centroid_s": None, "range_end_s": None, "is_on_bound": False}
    model_lines = [open(best_file).readline().split() for best_file in glob.glob(os.path.join(inversion_dir, "[0-9]*best.dat"))]
    snapshot = os.path.join(inversion_dir, "best_model_hybrid.dat")
    if not model_lines and os.path.exists(snapshot):
        model_lines = [open(snapshot).readline().split()]
    model_lines = [fields for fields in model_lines if len(fields) > 2]
    search_path = os.path.join(inversion_dir, "search_par.file")
    if not model_lines or not os.path.exists(search_path):
        return unknown
    range_bounds = {}
    for line in open(search_path):
        fields = line.split()
        if len(fields) >= 2 and fields[0] in ("cen_min", "cen_max"):
            range_bounds[fields[0]] = float(fields[1])
    if len(range_bounds) < 2:
        return unknown
    # a best-sample line: subevent count, "<neq", then 8 numbers per subevent: misfit, centroid time, ...
    best_fields = min(model_lines, key=lambda fields: float(fields[2]))
    last_centroid = max(float(best_fields[3 + 8 * subevent]) for subevent in range(int(best_fields[0])))
    room_left = range_bounds["cen_max"] - last_centroid
    return {"last_centroid_s": last_centroid, "range_end_s": range_bounds["cen_max"],
            "is_on_bound": room_left < bound_share_of_range * (range_bounds["cen_max"] - range_bounds["cen_min"])}


def island_best_misfits(inversion_dir):
    """best misfit of every island before finv's global resampling round (empty for v2 runs and for
    the folders of an example that ships without it)"""
    island_file = os.path.join(inversion_dir, "island_best_misfits.txt")
    if not os.path.exists(island_file):
        return np.array([])
    rows = np.loadtxt(island_file, ndmin=2)
    return np.sort(rows[:, 2]) if rows.size else np.array([])


def read_misfit_components(inversion_dir, best_misfit):
    """columns of misfit_components.txt (the factors of the best sample of every chain, written by
    evaluate_misfit_components.sh) as {name: array}, plus "worst_axis_deg" (the larger of the worst P and
    T rotations) and "best_index" (the row of the folder's best model). None when the file is missing,
    and also when its lowest sampled misfit is not the folder's best misfit: the file is then left over
    from an earlier run of the folder and describes other models."""
    components_path = os.path.join(inversion_dir, "misfit_components.txt")
    if not os.path.exists(components_path):
        return None
    table = np.atleast_1d(np.genfromtxt(components_path, names=True, dtype=None, encoding=None))
    columns = {name: np.asarray(table[name], dtype=float) for name in table.dtype.names if name != "chain"}
    if len(columns["sampled_misfit"]) == 0 or abs(columns["sampled_misfit"].min() / best_misfit - 1) > 1e-4:
        return None
    columns["worst_axis_deg"] = np.maximum(columns["worst_p_deg"], columns["worst_t_deg"])
    columns["best_index"] = int(np.argmin(columns["sampled_misfit"]))
    return columns


def akaike_criterion(subevent_counts, best_misfits):
    """AIC of every N, up to a constant: n ln(misfit) + 2 k N with n = 1 / sigma^2 (see the docstring)"""
    return np.log(best_misfits) / posterior_sigma**2 + 2 * parameters_per_subevent * np.asarray(subevent_counts)


# ---------------------------------------------------------------------------------------------------
# the best misfit of every N (all chains, end of the run), its AIC and the selection
# ---------------------------------------------------------------------------------------------------
rows = []
for inversion_dir in glob.glob("inv_*_*sub"):
    match = re.search(r"_(\d+)sub$", inversion_dir)
    if not match:
        continue
    chains = chain_best_misfits(inversion_dir)
    if len(chains) == 0:
        continue
    rows.append({"count": int(match.group(1)), "chains": chains, "islands": island_best_misfits(inversion_dir),
                 "components": read_misfit_components(inversion_dir, chains.min()),
                 "bound": last_centroid_against_its_bound(inversion_dir)})
rows.sort(key=lambda row: row["count"])
if len(rows) < 2:
    sys.exit("plot_lcurve: need results of at least two subevent counts (inv_*_<N>sub with *best.dat)")
counts = [row["count"] for row in rows]
best_per_n = np.array([row["chains"].min() for row in rows])
# an N whose best model ends on the bound of its centroid range cannot be selected (see the
# docstring); when every N does, nothing can be excluded
is_on_bound = np.array([row["bound"]["is_on_bound"] for row in rows])
every_count_is_on_bound = bool(is_on_bound.all())
can_be_selected = np.ones(len(rows), dtype=bool) if every_count_is_on_bound else ~is_on_bound
aic_above_lowest = akaike_criterion(counts, best_per_n)
aic_above_lowest -= aic_above_lowest[can_be_selected].min()
lowest_aic_count = min((difference, count) for difference, count, allowed in zip(aic_above_lowest, counts, can_be_selected) if allowed)[1]
# every N that the data cannot tell apart from the one with the lowest AIC (that one included); the
# selection is the smallest of them
tied_counts = [count for count, difference, allowed in zip(counts, aic_above_lowest, can_be_selected)
               if allowed and difference < tie_aic_difference]
selected = tied_counts[0]
selected_index = counts.index(selected)

# one row per N, every value right-aligned under its column name (the parsers of the tools read the
# first two columns as "a number, then the best misfit", whatever the spacing)
column_names = ["n", "best_misfit", "reduction_%", "delta_aic", "chains", "chains_within_1%"]
report = ["  ".join(column_names) + "  island_bests"]
for index, row in enumerate(rows):
    # how much lower than the N before it (negative: the misfit rose, the search found no place where
    # one more subevent helps)
    reduction_text = "-" if index == 0 else f"{100 * (1 - best_per_n[index] / best_per_n[index - 1]):.1f}"
    near_best = int(np.sum(row["chains"] <= 1.01 * best_per_n[index]))
    island_text = ",".join(f"{value:.4f}" for value in row["islands"][:8]) or "-"
    aic_text = f"{aic_above_lowest[index]:.1f}" if can_be_selected[index] else "bound"
    values = [str(row["count"]), f"{best_per_n[index]:.4f}", reduction_text, aic_text, str(len(row["chains"])), str(near_best)]
    report.append("  ".join(value.rjust(len(name)) for value, name in zip(values, column_names)) + "  " + island_text)
# the factors of each best model; rows start with "N=", which the parsers of the table above
# (tools/compare_with_reference.py, tools/make_readme_figures.py: a number, then the best misfit) skip
component_rows = [row for row in rows if row["components"] is not None]
if component_rows:
    report.append("misfit components of each best model: penalized = data_misfit x clvd_factor x stress_factor "
                  "(evaluate_misfit_components.sh)")
    component_names = ["penalized", "data_misfit", "clvd_factor", "stress_factor", "worst_axis_deg", "p_misfit",
                       "sh_misfit", "regional_misfit"]
    report.append("n    " + "  ".join(component_names))
    for row in component_rows:
        components, best_index = row["components"], row["components"]["best_index"]
        value = {name: components[name][best_index] for name in components if name != "best_index"}
        values = [f"{value['penalized_misfit']:.4f}", f"{value['data_misfit']:.4f}", f"{value['clvd_factor']:.4f}",
                  f"{value['stress_factor']:.4f}", f"{value['worst_axis_deg']:.1f}", f"{value['p_misfit']:.3f}",
                  f"{value['sh_misfit']:.3f}", f"{value['regional_misfit']:.3f}"]
        report.append(f"N={row['count']:<3d}" + "  ".join(text.rjust(len(name)) for text, name in zip(values, component_names)))
counts_without_components = [str(row["count"]) for row in rows if row["components"] is None]
if counts_without_components:
    report.append(f"NOTE: no misfit components for N={','.join(counts_without_components)} (misfit_components.txt "
                  f"missing or left over from an earlier run; step 3 writes it since 2026-11-16)")
report.append(f"selection: AIC = ln(misfit) / sigma^2 + 2 k N, sigma = {posterior_sigma}, k = {parameters_per_subevent} "
              f"parameters per subevent, so a subevent costs {100 * log_misfit_drop_per_subevent:g}% of misfit (ln(misfit) "
              f"has to fall by {log_misfit_drop_per_subevent:g}, the misfit by {100 * (1 - np.exp(-log_misfit_drop_per_subevent)):.1f}%, "
              f"to tie); delta_aic = AIC above the lowest; every N with delta_aic below {tie_aic_difference:g} is tied, "
              f"and the smallest of them is selected")
if len(tied_counts) > 1:
    tie_text = ", ".join(f"N={count} (delta_aic {aic_above_lowest[counts.index(count)]:.1f})" for count in tied_counts)
    report.append(f"NOTE: a tie - {tie_text} are equally supported (delta_aic below {tie_aic_difference:g}); "
                  f"N={selected} is selected as the smallest: report {'both' if len(tied_counts) == 2 else 'all of them'}")
bound_rows = [row for row in rows if row["bound"]["is_on_bound"]]
if bound_rows:
    bound_text = (f"the last subevent of the best model lies on the upper bound of the centroid range (N={bound_rows[0]['count']}: "
                  f"{bound_rows[0]['bound']['last_centroid_s']:.1f} of {bound_rows[0]['bound']['range_end_s']:g} s), where it fits later moment release, or "
                  f"arrivals that the model does not predict, from the edge of the range. If the rupture may have lasted longer, "
                  f"raise cen_max and the windows and run step 3 again (README sections 4.6 and 11)")
    if every_count_is_on_bound:
        report.append(f"NOTE: for every N {bound_text}. No N could be excluded: do not trust the selection below")
    else:
        report.append(f"NOTE: N={','.join(str(row['count']) for row in bound_rows)} cannot be selected (\"bound\" in the column "
                      f"delta_aic): {bound_text}")
if lowest_aic_count == counts[-1]:
    # ./submit extends the grid while the lowest AIC sits at its edge, up to its limit of 7 subevents;
    # an AIC that is still lowest at the edge is where the search stopped, not an optimum (v3 Kamchatka
    # rerun: the misfit fell another 10% from 7 to 8 subevents)
    report.append(f"NOTE: the AIC is lowest at N={lowest_aic_count}, the largest number of subevents tried; more "
                  f"subevents might fit better, so read {selected} as a lower bound")
report.append(f"lowest AIC: N={lowest_aic_count}")
report.append(f"SUGGESTED number of subevents (smallest N with delta_aic below {tie_aic_difference:g}): {selected}")
open("lcurve.txt", "w").write("\n".join(report) + "\n")
print("\n".join(report))

# ---------------------------------------------------------------------------------------------------
# figure: the L-curve (left) and, when step 3 evaluated them, its factors (right, 2 x 2). One search,
# one colour. Left: line = best misfit per N, star = the selection, dashed line = the misfit with
# the AIC of the selection (5% lower per added subevent): the L-curve touches it at the star and
# lies above it elsewhere, except that an N tied with the star can lie up to 0.5% below it (2 units
# of AIC), and so can an N that cannot be selected (hollow point: its last subevent lies on the bound
# of the centroid range); small open circles = the island bests. Right: line = the factors of the
# same best models, dots = the best sample of each chain.
# ---------------------------------------------------------------------------------------------------
line_colour = "#2a78d6"
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6})
if component_rows:
    figure = plt.figure(figsize=(9, 4.4))
    panel_grid = figure.add_gridspec(2, 3, width_ratios=[1.25, 1, 1], wspace=0.42, hspace=0.12)
    lcurve_axis = figure.add_subplot(panel_grid[:, 0])
    factor_axes = {key: figure.add_subplot(panel_grid[grid_row, grid_column], sharex=lcurve_axis)
                   for key, (grid_row, grid_column) in zip(factor_labels, ((0, 1), (0, 2), (1, 1), (1, 2)))}
else:
    figure, lcurve_axis = plt.subplots(figsize=(6, 4.5))
    factor_axes = {}

selected_misfit = best_per_n[selected_index]
tie_line_counts = np.linspace(counts[0] - 0.3, counts[-1] + 0.3, 50)
tie_line_misfits = selected_misfit * np.exp(-log_misfit_drop_per_subevent * (tie_line_counts - selected))
lcurve_axis.plot(tie_line_counts, tie_line_misfits, "--", color="0.55", lw=0.8, zorder=1)
for row in rows:
    lcurve_axis.plot([row["count"]] * len(row["islands"]), row["islands"], "o", ms=3.5, mfc="none",
                     mec=line_colour, mew=0.7, alpha=0.6)
lcurve_axis.plot(counts, best_per_n, "-o", color=line_colour, lw=1.4, ms=4.5)
if bound_rows and not every_count_is_on_bound:
    # N that cannot be selected (their last subevent lies on the bound of the centroid range): hollow points
    bound_counts = [row["count"] for row in bound_rows]
    bound_misfits = [best_per_n[counts.index(count)] for count in bound_counts]
    lcurve_axis.plot(bound_counts, bound_misfits, "o", ms=5.5, mfc="white", mec=line_colour, mew=1.2, zorder=4)
lcurve_axis.plot(selected, selected_misfit, "*", ms=14, color=line_colour, mec="white", mew=0.8, zorder=5)
# island bests far above the curve (two California 6-subevent islands stuck at 5.3 and 5.7) are left
# outside the axes, where they would flatten the curve itself
lcurve_axis.set_ylim(min(best_per_n.min(), tie_line_misfits.min()) * 0.97, best_per_n.max() * 1.05)
lcurve_axis.set_ylabel("penalized misfit")

random_generator_for_jitter = np.random.default_rng(1)       # fixed: the chain dots do not move between runs
for key, axis in factor_axes.items():
    for row in component_rows:
        chain_values = row["components"][key]
        jitter = random_generator_for_jitter.uniform(-0.08, 0.08, len(chain_values))
        axis.scatter(row["count"] + jitter, chain_values, s=5, color=line_colour, alpha=0.3, lw=0)
    axis.plot([row["count"] for row in component_rows],
              [row["components"][key][row["components"]["best_index"]] for row in component_rows],
              "-o", color=line_colour, lw=1.4, ms=4.5)
    axis.set_ylabel(factor_labels[key])
if factor_axes:
    factor_axes["clvd_factor"].axhline(1.0, color="0.6", lw=0.6, zorder=0)
    factor_axes["stress_factor"].axhline(1.0, color="0.6", lw=0.6, zorder=0)
    largest_stress_factor = max(np.nanmax(row["components"]["stress_factor"]) for row in component_rows)
    factor_axes["stress_factor"].set_ylim(0.995, max(1.05, 1.01 * largest_stress_factor))
    factor_axes["worst_axis_deg"].axhline(stress_cone_deg, color="0.6", lw=0.8, ls="--", zorder=0)
    largest_axis_rotation = max(np.nanmax(row["components"]["worst_axis_deg"]) for row in component_rows)
    factor_axes["worst_axis_deg"].set_ylim(0, max(40, largest_axis_rotation + 5))
    for key in ("data_misfit", "clvd_factor"):
        plt.setp(factor_axes[key].get_xticklabels(), visible=False)
    for key in ("stress_factor", "worst_axis_deg"):
        factor_axes[key].set_xlabel("subevents")
for axis in [lcurve_axis, *factor_axes.values()]:
    axis.grid(axis="y", ls="--", lw=0.4, color="0.85")
    axis.set_axisbelow(True)
lcurve_axis.set_xticks(counts)
lcurve_axis.set_xlim(counts[0] - 0.3, counts[-1] + 0.3)
lcurve_axis.set_xlabel("subevents")
# the label of the dashed line runs along it; its angle on the page is known only now that the axis
# limits are final. It is written 2 points below the line, where the L-curve never is (above it, the
# label sat on the first points of the Vanuatu curve), on the side of the star that has more room:
# from the left end of the line, or just right of the star when the selection is among the first N
line_start_on_page, line_end_on_page = lcurve_axis.transData.transform(
    [(tie_line_counts[0], tie_line_misfits[0]), (tie_line_counts[-1], tie_line_misfits[-1])])
line_angle = np.arctan2(*(line_end_on_page - line_start_on_page)[::-1])
star_is_nearer_the_left_end = selected - counts[0] < counts[-1] - selected
label_count = selected + 0.3 if star_is_nearer_the_left_end else tie_line_counts[2]
label_misfit = selected_misfit * np.exp(-log_misfit_drop_per_subevent * (label_count - selected))
lcurve_axis.annotate(f"{100 * log_misfit_drop_per_subevent:g}% per subevent", xy=(label_count, label_misfit),
                     xytext=(2 * np.sin(line_angle), -2 * np.cos(line_angle)), textcoords="offset points",
                     rotation=np.degrees(line_angle), rotation_mode="anchor", ha="left", va="top", fontsize=7, color="0.4")
figure.savefig("lcurve.pdf", bbox_inches="tight")
print("wrote lcurve.pdf / lcurve.txt")
