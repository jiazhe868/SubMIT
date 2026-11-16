#!/usr/bin/env python3
"""L-curve (best penalized misfit against the number of subevents), the selected number, and how
reproducible that selection is.

Rule (README section 4.6, unchanged since 2026-08-11): the SMALLEST N whose best misfit is within 5%
of the lowest misfit over the whole grid. Scale-free and anchored to the global minimum.

Why a reproducibility estimate. The best misfit of each N comes out of a stochastic search, so the
selection is only as reliable as the search. In the 2026 Panama clean-environment run (v2 code) the
selected N = 2 rested on ONE of 16 chains (1.530): the median chain's best was 17% worse, and
without that chain (next best 1.693) the rule would have picked 3. The best misfits of independent
replicate searches are therefore bootstrapped (resampled with replacement, 1000 times, fixed seed)
and the share of resamples that select each N is reported; a winner below 80% is flagged as
unstable.

What counts as a replicate. v3 chains are not independent: finv resamples them within islands of 8
chains and, in the last burn-in round, across all chains, after which every chain descends from the
best island (finv.f90). Their bests then agree within 1-2% even when two seeds of the same run differ
by 5% (California 3sub: 1.761 and 1.847), so a bootstrap over chains would look stable when the search
is not. The replicates are therefore the islands, whose bests before the global round finv writes to
island_best_misfits.txt; with 32 ranks there are 4 islands, with 16 ranks 2, and with 8 or fewer one,
which leaves the reproducibility unmeasured. An island best is the lowest misfit the island visited
during the burn-in; best_misfit is the best post-burn-in sample of all chains, the model that is
reported. An island best can therefore lie slightly below best_misfit: by up to 0.5% in the v3 reruns
of the four examples (Venezuela 6sub: 2.0877 against 2.0976). Folders without that file (v2 runs, independent chains
with anneal_stages 1) fall back to the chains themselves. The published
examples ship chains for the selected N only, so their other N count as one replicate, and a selection
that rests on single replicates is reported as not measured.

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

selection_band = 0.05
bootstrap_count = 1000
stable_share = 0.8


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


def island_best_misfits(inversion_dir):
    """best misfit of every island before finv's global resampling round (empty if not written)"""
    island_file = os.path.join(inversion_dir, "island_best_misfits.txt")
    if not os.path.exists(island_file):
        return np.array([])
    rows = np.loadtxt(island_file, ndmin=2)
    return rows[:, 2] if rows.size else np.array([])


def select(best_per_n, counts):
    """smallest N within the band of the grid minimum"""
    grid_minimum = min(best_per_n)
    for count, best in zip(counts, best_per_n):
        if best <= (1 + selection_band) * grid_minimum:
            return count
    return counts[-1]


# ---------------------------------------------------------------------------------------------------
# results of every N: the overall best (all chains, end of the run) and the replicate bests
# ---------------------------------------------------------------------------------------------------
rows = []
for inversion_dir in glob.glob("inv_*_*sub"):
    match = re.search(r"_(\d+)sub$", inversion_dir)
    if not match:
        continue
    chains = chain_best_misfits(inversion_dir)
    if len(chains) == 0:
        continue
    if os.path.exists(os.path.join(inversion_dir, "island_best_misfits.txt")):
        # a v3 run: its chains are copies of the islands' descendants, so only islands count, even
        # when there is a single one (8 or fewer ranks: reproducibility then shows as not measured)
        islands = island_best_misfits(inversion_dir)
        replicates, replicate_kind = (islands if len(islands) else np.array([chains.min()])), "islands"
    else:
        replicates, replicate_kind = chains, "chains"
    rows.append({"count": int(match.group(1)), "chains": chains, "replicates": replicates, "kind": replicate_kind})
rows.sort(key=lambda row: row["count"])
if len(rows) < 2:
    sys.exit("plot_lcurve: need results of at least two subevent counts (inv_*_<N>sub with *best.dat)")
counts = [row["count"] for row in rows]
best_per_n = [row["chains"].min() for row in rows]
selected = select(best_per_n, counts)

# ---------------------------------------------------------------------------------------------------
# bootstrap over replicates: how often would a search with one replicate's effort pick each N?
# ---------------------------------------------------------------------------------------------------
random_generator = np.random.default_rng(12345)       # fixed seed: the report is reproducible
votes = {count: 0 for count in counts}
for _ in range(bootstrap_count):
    resampled_best = [random_generator.choice(row["replicates"], size=len(row["replicates"]), replace=True).min()
                      for row in rows]
    votes[select(resampled_best, counts)] += 1
share = votes[selected] / bootstrap_count
# every N enters the rule (the smaller N through their distance from the band), so every N needs
# replicates; the published examples (chains of the selected N only) are therefore not measured
is_measured = all(len(row["replicates"]) >= 2 for row in rows)

report = ["n  best_misfit  chains  chains_within_1%  replicates  replicate_bests  share_selected_in_bootstrap"]
for row, best in zip(rows, best_per_n):
    near_best = int(np.sum(row["chains"] <= 1.01 * best))
    replicate_text = ",".join(f"{value:.4f}" for value in np.sort(row["replicates"])[:8])
    report.append(f"{row['count']}  {best:.4f}  {len(row['chains'])}  {near_best}  {len(row['replicates'])} {row['kind']}  "
                  f"{replicate_text}  {votes[row['count']] / bootstrap_count:.2f}")
for index in range(1, len(counts)):
    if best_per_n[index] > best_per_n[index - 1] + 1e-4:
        report.append(f"NOTE: penalized misfit rises {counts[index - 1]}sub->{counts[index]}sub "
                      f"({best_per_n[index - 1]:.4f}->{best_per_n[index]:.4f}) - expected when the stress/CLVD "
                      f"penalties reject extra subevents")
grid_minimum = min(best_per_n)
report.append(f"grid minimum {grid_minimum:.4f}; 5% band <= {(1 + selection_band) * grid_minimum:.4f}")
if not is_measured:
    report.append("NOTE: reproducibility not measured - the deciding N have a single replicate each")
elif share < stable_share:
    runner_up = max((count for count in counts if count != selected), key=lambda count: votes[count])
    report.append(f"WARNING: unstable selection - N={selected} wins only {100 * share:.0f}% of the replicate bootstrap "
                  f"(N={runner_up}: {100 * votes[runner_up] / bootstrap_count:.0f}%). The deciding best misfits differ "
                  f"between independent searches: rerun those N with more MPI ranks (SUBMIT_NCORE) or another "
                  f"rank.dat offset, and report both candidates.")
if selected == counts[-1]:
    # ./submit extends the grid while the selection sits at its edge, up to its limit of 7 subevents;
    # a selection still at the edge is where the search stopped, not an elbow (v3 Kamchatka rerun:
    # the misfit fell another 10% from 7 to 8 subevents)
    report.append(f"NOTE: N={selected} is the largest number of subevents tried; more subevents might fit "
                  f"better, so read {selected} as a lower bound")
report.append(f"SUGGESTED number of subevents (elbow): {selected}")
open("lcurve.txt", "w").write("\n".join(report) + "\n")
print("\n".join(report))

# ---------------------------------------------------------------------------------------------------
# figure: best misfit per N, replicate bests as dots
# ---------------------------------------------------------------------------------------------------
figure, axis = plt.subplots(figsize=(6, 4.5))
axis.plot(counts, best_per_n, "o-", color="#1f4e8c")
for row in rows:
    axis.plot([row["count"]] * len(row["replicates"]), row["replicates"], ".", color="#1f4e8c", alpha=0.35, markersize=5)
share_text = f"{100 * share:.0f}% of bootstrap" if is_measured else "reproducibility not measured"
axis.axvline(selected, color="crimson", ls="--", lw=1, label=f"selected: {selected} ({share_text})")
axis.set_ylim(top=max(best_per_n) * 1.3)
axis.set_xlabel("number of subevents")
axis.set_ylabel("best misfit")
axis.set_xticks(counts)
axis.grid(alpha=0.3)
axis.legend()
figure.tight_layout()
figure.savefig("lcurve.pdf")
print("wrote lcurve.pdf / lcurve.txt")
