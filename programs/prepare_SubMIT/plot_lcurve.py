#!/usr/bin/env python3
"""L-curve (best misfit vs number of subevents) with monotonicity checks and an
elbow suggestion. Run in IRIS/ after step3. Writes lcurve.pdf + lcurve.txt.
Uses best_model_exploration.dat snapshots when present (reseeded ensemble reruns
overwrite *best.dat, so raw scans can pick up herded ensemble values)."""
import glob, os, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = []
for d in sorted(glob.glob("inv_*_?sub")):
    n = int(re.search(r"_(\d)sub$", d).group(1))
    # pick the NEWEST best_model_* snapshot (stale snapshots from earlier
    # runs under different weights are NOT comparable - 2026-08-11 Chile:
    # old exploration files silently mixed into a hybrid-mode curve)
    cands = [(os.path.getmtime(f), f, m) for m in
             ("hybrid", "ensemble", "exploration")
             for f in [os.path.join(d, f"best_model_{m}.dat")]
             if os.path.exists(f)]
    if cands:
        _, snap, src = max(cands)
        src = f"{src}-snapshot"
        best = float(open(snap).read().split()[2])
    else:
        vals = []
        for f in glob.glob(os.path.join(d, "*best.dat")):
            t = open(f).readline().split()
            if len(t) > 2:
                vals.append(float(t[2]))
        if not vals:
            continue
        best = min(vals)
        src = "raw *best.dat (no exploration snapshot - may be an ensemble run)"
    rows.append((n, best, src))
rows.sort()
if len(rows) < 2:
    sys.exit("plot_lcurve: need >=2 inv_*_Nsub results")
ns = np.array([r[0] for r in rows]); ms = np.array([r[1] for r in rows])

report = ["n  best_misfit  source"]
report += [f"{n}  {m:.4f}  {s}" for n, m, s in rows]
for i in range(1, len(ns)):
    if ms[i] > ms[i-1] + 1e-4:
        report.append(f"NOTE: penalized misfit rises {ns[i-1]}sub->{ns[i]}sub "
                      f"({ms[i-1]:.4f}->{ms[i]:.4f}) - expected when the stress/CLVD "
                      f"penalties reject extra subevents; suspect non-convergence only "
                      f"if the DATA misfit also rises")
# elbow (user rule 2026-08-11): the SMALLEST n whose best misfit is within
# 5% of the lowest misfit found across the whole subevent-count grid -
# scale-free and anchored to the global minimum, unlike slope heuristics
m_min = min(ms)
elbow = ns[-1]
for i in range(len(ns)):
    if ms[i] <= 1.05 * m_min:
        elbow = ns[i]
        break
report.append(f"grid minimum {m_min:.4f}; 5% band <= {1.05*m_min:.4f}")
report.append(f"SUGGESTED number of subevents (elbow): {elbow}")
open("lcurve.txt", "w").write("\n".join(report) + "\n")
print("\n".join(report))

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(ns, ms, "o-", color="#1f4e8c")
ax.axvline(elbow, color="crimson", ls="--", lw=1, label=f"elbow: {elbow} subevents")
ax.set_xlabel("number of subevents"); ax.set_ylabel("best misfit")
ax.set_xticks(ns); ax.grid(alpha=0.3); ax.legend()
fig.tight_layout(); fig.savefig("lcurve.pdf")
print("wrote lcurve.pdf / lcurve.txt")
