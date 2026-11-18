#!/usr/bin/env python3
"""Pool the Markov chains of one inversion folder into a posterior ensemble and report convergence.

What it writes (run inside an inv_<event>_<N>sub folder after finv):
  pooled_ensemble.dat    the post-burn-in samples of the converged chains (finv's chain format)
  allsamples.dat         misfit, then (centroid time, x, y, duration, depth) per subevent, one
                         sample per row; the 95% error bars of the subevent map are read from it
  histoplot_py.pdf       posterior histograms, one row per subevent, with 95% intervals
  misfit_evolution.pdf/.png   misfit along every chain (convergence check)
  convergence.txt        chains near the best, split R-hat of every parameter (README section 9)

Why a convergence gate. A chain whose best misfit stays far above the global best never reached
the main posterior mode: with sigma = 0.05 a chain 5% above the best sits at exp(-10) of the
density. Such chains are left out (gate: chain best <= 1.05 x global best). The included chains
are exchangeable draws of one posterior and are pooled unweighted: weighting them by misfit would
count the likelihood twice.

History: replaces pool_chains.py, process_and_plot.py and plot_misfit_evolution.py (v2). The v2
histogram script ranked chains by the LAST line of each *best.dat, which finv writes best-first,
i.e. by each chain's 50th-best sample.

usage: python posterior_summary.py [gate]      (default gate 0.05)
"""
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# reading finv's chain files: "NN <neq " then 8 numbers per subevent in 10-character fields:
# misfit, centroid time, x east, y north, duration, rupture velocity, rupture direction, depth
# ---------------------------------------------------------------------------


def read_chain_file(path):
    """rows of floats, one per sample; fixed-width parsing survives fields that touch"""
    rows = []
    for line in open(path):
        if "<neq" not in line:
            continue
        body = line.split("<neq", 1)[1].rstrip("\n")
        body = body[1:] if body.startswith(" ") else body
        fields = [body[i:i + 10] for i in range(0, len(body) - 9, 10)]
        try:
            rows.append([float(field) for field in fields])
        except ValueError:
            continue
    return np.array(rows)


def split_r_hat(chains):
    """Gelman-Rubin statistic on half-chains (split R-hat, Vehtari et al. 2021, without rank
    normalization). Values near 1 mean the chains sample the same distribution; above ~1.1 they
    disagree. Returns nan with fewer than two usable chains."""
    halves = []
    for chain in chains:
        half_length = len(chain) // 2
        if half_length >= 10:
            halves += [chain[:half_length], chain[half_length:2 * half_length]]
    if len(halves) < 4:
        return float("nan")
    length = min(len(half) for half in halves)
    stacked = np.array([half[:length] for half in halves])
    within = stacked.var(axis=1, ddof=1).mean()
    between = length * stacked.mean(axis=1).var(ddof=1)
    if within <= 0:
        return float("nan")
    return float(np.sqrt(((length - 1) / length * within + between / length) / within))


# ---------------------------------------------------------------------------
# gate and pool
# ---------------------------------------------------------------------------
gate_relative = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
chain_files = sorted(glob.glob("[0-9]*chain.dat"), key=lambda name: int(re.match(r"(\d+)", name).group(1)))
chains = {}
for chain_file in chain_files:
    rows = read_chain_file(chain_file)
    if len(rows):
        chains[chain_file] = rows
if not chains:
    sys.exit("posterior_summary: no <rank>chain.dat with samples here")
subevent_count = next(iter(chains.values())).shape[1] // 8
chain_best = {name: rows[:, 0].min() for name, rows in chains.items()}
global_best = min(chain_best.values())
kept_chains = [name for name in chains if chain_best[name] <= (1 + gate_relative) * global_best]
with open("pooled_ensemble.dat", "w") as pooled:
    for name in kept_chains:
        pooled.writelines(open(name).readlines())
print(f"posterior_summary: global best {global_best:.4f}; gate <= {(1 + gate_relative) * global_best:.4f}; "
      f"kept {len(kept_chains)}/{len(chains)} chains")

# samples: misfit, then per subevent (centroid time, x, y, duration, depth); subevents in time order
pooled_rows = np.vstack([chains[name] for name in kept_chains])
columns = [0]
for subevent in range(subevent_count):
    columns += [8 * subevent + 1, 8 * subevent + 2, 8 * subevent + 3, 8 * subevent + 4, 8 * subevent + 7]
samples = pooled_rows[:, columns]
for row in samples:
    blocks = row[1:].reshape(subevent_count, 5)
    row[1:] = blocks[np.argsort(blocks[:, 0])].ravel()
np.savetxt("allsamples.dat", samples, fmt="%.6f")

# ---------------------------------------------------------------------------
# convergence report
# ---------------------------------------------------------------------------
parameter_names = ["centroid time", "x east", "y north", "duration", "depth"]
report = [f"chains: {len(chains)}; best misfit {global_best:.4f}"]
for tolerance in (0.01, 0.02, 0.05):
    within = sum(1 for best in chain_best.values() if best <= (1 + tolerance) * global_best)
    report.append(f"chains whose best is within {100 * tolerance:.0f}% of the global best: {within}")
report.append("split R-hat of the gated chains (1 = converged, > 1.1 = chains disagree):")
largest_r_hat = 0.0
for subevent in range(subevent_count):
    values = []
    for parameter, column in zip(parameter_names, (1, 2, 3, 4, 7)):
        values_of_chains = [chains[name][:, 8 * subevent + column] for name in kept_chains]
        if np.ptp(np.concatenate(values_of_chains)) < 1e-6:
            continue                       # not sampled: E1 sits at the hypocentre, epicentre and depth (finv.f90)
        r_hat = split_r_hat(values_of_chains)
        largest_r_hat = max(largest_r_hat, r_hat) if np.isfinite(r_hat) else largest_r_hat
        values.append(f"{parameter} {r_hat:.2f}")
    report.append(f"  E{subevent + 1}: " + ", ".join(values))
report.append(f"largest split R-hat: {largest_r_hat:.2f}")
open("convergence.txt", "w").write("\n".join(report) + "\n")
print("\n".join(report))

# ---------------------------------------------------------------------------
# histograms: one row per subevent, shared range per column, 95% interval drawn at the base
# ---------------------------------------------------------------------------
plt.rcParams.update({"font.size": 12})
figure, axes = plt.subplots(max(subevent_count, 2), 5, figsize=(15, 2.8 * max(subevent_count, 2)))
# panels in the published layout; numbers index the (time, x, y, duration, depth) block of a subevent
panel_order = [0, 3, 1, 2, 4]
panel_titles = ["Centroid time (s)", "Duration (s)", "WE location (km)", "NS location (km)", "Depth (km)"]
for subevent in range(subevent_count):
    for panel, parameter_index in enumerate(panel_order):
        axis = axes[subevent, panel]
        values = samples[:, 1 + 5 * subevent + parameter_index]
        all_subevents = samples[:, [1 + 5 * other + parameter_index for other in range(subevent_count)]].ravel()
        padding = 0.06 * (all_subevents.max() - all_subevents.min() + 1e-9)
        axis_range = (all_subevents.min() - padding, all_subevents.max() + padding)
        if values.max() - values.min() < 1e-3:
            # a pinned parameter (E1 at the epicentre) would be a single 1e6-high bin
            axis.axvline(values[0], color="steelblue", lw=2)
            axis.annotate(f"fixed at {values[0]:.2f}", (0.5, 0.5), xycoords="axes fraction", ha="center", fontsize=10)
        else:
            # ~25 bins over the marginal's own span, snapped to 1/2/2.5/5 x 10^k (arbitrary widths beat
            # against the sampler's discreteness), never finer than 1/120 of the shared axis
            raw_width = max((values.max() - values.min()) / 25.0, (axis_range[1] - axis_range[0]) / 120.0)
            magnitude = 10.0 ** np.floor(np.log10(raw_width))
            bin_width = min((factor * magnitude for factor in (1.0, 2.0, 2.5, 5.0, 10.0)),
                            key=lambda width: abs(width - raw_width) if width >= raw_width else 1e9)
            edges = np.arange(np.floor(values.min() / bin_width) * bin_width,
                              np.ceil(values.max() / bin_width) * bin_width + 1.5 * bin_width, bin_width)
            edge_color = "black" if bin_width >= (axis_range[1] - axis_range[0]) / 60.0 else "none"
            axis.hist(values, bins=edges, density=True, edgecolor=edge_color, facecolor="skyblue", alpha=0.85)
            lower, upper = np.percentile(values, [2.5, 97.5])
            axis.plot([lower, upper], [0, 0], "-ks", linewidth=2, markersize=5, markerfacecolor="gray")
        axis.set_xlim(axis_range)
        if panel == 0:
            axis.set_ylabel(f"E{subevent + 1} prob.")
        if subevent == subevent_count - 1:
            axis.set_xlabel(panel_titles[panel])
plt.tight_layout()
plt.savefig("histoplot_py.pdf", dpi=300)
plt.close()

# ---------------------------------------------------------------------------
# misfit along every chain, with the running ensemble best
# ---------------------------------------------------------------------------
figure, axis = plt.subplots(figsize=(8.5, 5.5))
colors = plt.cm.viridis(np.linspace(0, 0.9, len(chains)))
longest = max(len(rows) for rows in chains.values())
padded = np.full((len(chains), longest), np.nan)
for index, (color, (name, rows)) in enumerate(zip(colors, chains.items())):
    axis.plot(np.arange(1, len(rows) + 1), rows[:, 0], lw=0.6, alpha=0.55, color=color)
    padded[index, :len(rows)] = rows[:, 0]
    padded[index, len(rows):] = rows[-1, 0]
axis.plot(np.arange(1, longest + 1), np.minimum.accumulate(np.nanmin(padded, axis=0)), "k", lw=1.8,
          label="ensemble best (running)")
axis.set_xlabel("sample number (post burn-in)", fontsize=12)
axis.set_ylabel("misfit", fontsize=12)
axis.set_title(f"Misfit evolution, {len(chains)} chains", fontsize=12)
axis.legend(fontsize=10, frameon=False)
axis.grid(alpha=0.25)
figure.tight_layout()
figure.savefig("misfit_evolution.pdf")
figure.savefig("misfit_evolution.png", dpi=140)
print(f"posterior_summary: wrote allsamples.dat ({len(samples)} samples), histoplot_py.pdf, "
      f"misfit_evolution.pdf/.png, convergence.txt")
