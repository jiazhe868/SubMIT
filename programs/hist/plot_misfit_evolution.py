#!/usr/bin/env python3
"""MCMC convergence diagnostic: misfit vs sample number for every chain.

Run inside an inv_* directory after finv. Reads <rank>chain.dat (post-burn-in
samples; token 2 of each row is the misfit) and overlays all chains as thin
lines plus the running ensemble minimum. Writes misfit_evolution.pdf/.png.
"""
import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

chains = sorted(glob.glob("[0-9]*chain.dat"),
                key=lambda f: int(re.match(r"(\d+)", f).group(1)))
if not chains:
    raise SystemExit("plot_misfit_evolution: no <rank>chain.dat files found")

series = []
for f in chains:
    mis = [float(line.split()[2]) for line in open(f) if line.split()]
    if mis:
        series.append((f, np.asarray(mis)))

n = max(len(m) for _, m in series)
fig, ax = plt.subplots(figsize=(8.5, 5.5))
cmap = plt.cm.viridis(np.linspace(0, 0.9, len(series)))
for c, (f, mis) in zip(cmap, series):
    ax.plot(np.arange(1, len(mis) + 1), mis, lw=0.6, alpha=0.55, color=c)
# running ensemble minimum across chains (pad shorter chains with their last value)
stack = np.full((len(series), n), np.nan)
for i, (_, mis) in enumerate(series):
    stack[i, :len(mis)] = mis
    stack[i, len(mis):] = mis[-1]
ax.plot(np.arange(1, n + 1), np.minimum.accumulate(np.nanmin(stack, axis=0)),
        "k", lw=1.8, label="ensemble best (running)")
ax.set_xlabel("sample number (post burn-in)", fontsize=12)
ax.set_ylabel("misfit", fontsize=12)
ax.set_title(f"Misfit evolution, {len(series)} chains", fontsize=12)
ax.legend(fontsize=10, frameon=False)
ax.grid(alpha=0.25)
ax.tick_params(labelsize=11)
fig.tight_layout()
fig.savefig("misfit_evolution.pdf")
fig.savefig("misfit_evolution.png", dpi=140)
best = np.nanmin(stack)
print(f"wrote misfit_evolution.pdf/.png ({len(series)} chains, "
      f"{n} samples/chain, best {best:.4f})")
