#!/usr/bin/env python3
"""Write figs_and_results/RESULTS.md: formatted subevent-inversion results for
the selected subevent count. usage: make_results_doc.py <event> <nsub> (in IRIS/)"""
import glob, os, sys
import numpy as np

ev, nn = sys.argv[1], int(sys.argv[2])
inv, fwd = f"inv_{ev}_{nn}sub", f"fwd_{ev}_{nn}sub"
model = np.loadtxt(f"{fwd}/Input.model", ndmin=2)
fm = np.loadtxt(f"{fwd}/fm.dat", ndmin=2)[:, 1:7]
def mw(m): return 2/3*np.log10(np.sqrt((m[0]**2+2*m[1]**2+2*m[2]**2+m[3]**2+2*m[4]**2+m[5]**2)/2)*1e27)-10.7
misfit = open(f"{fwd}/misfit.dat").read().strip()
mode = "unknown"
for m in ("hybrid", "ensemble", "exploration"):
    if os.path.exists(f"{inv}/best_model_{m}.dat"):
        mode = m; break
mwcat = ""
try:
    mwcat = f" (catalog Mw {float(open(f'{inv}/mainshock.dat').read().split()[3]):.1f})"
except Exception:
    pass
L = [f"# Subevent inversion results: {ev}", "",
     f"- Selected subevents: **{nn}** | best misfit **{misfit}** | sampler mode: {mode}",
     f"- Total moment: **Mw {mw(fm.sum(axis=0)):.2f}**{mwcat}", ""]
if os.path.exists("lcurve.txt"):
    L += ["## L-curve (misfit vs subevents)", "```", open("lcurve.txt").read().rstrip(), "```", ""]
L += ["## Best model", "",
      "| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw |",
      "|---|---|---|---|---|---|---|"]
for i, r in enumerate(model):
    L.append(f"| E{i+1} | {r[0]:.2f} | {r[1]:.2f} | {r[2]:.2f} | {r[3]:.2f} | {r[6]:.1f} | {mw(fm[i]):.2f} |")
L += ["", "## Figures", "",
      "- `subevents.pdf/png` — map (lon/lat, 95% error bars), depth section, STFs",
      "- `fits_P.pdf`, `fits_Pvel.pdf`, `fits_SH.pdf`, `fits_rayl.pdf` — waveform fits",
      "- `misfit_evolution.pdf` — misfit vs sample number, all chains",
      "- `histoplot_py.pdf` — posterior histograms (gated ensemble)",
      "- `lcurve.pdf` — model selection", ""]
try:
    b = [float(open(f).readline().split()[2]) for f in glob.glob(f"{inv}/*best.dat")]
    g = min(b); w = sum(1 for x in b if x <= 1.1*g)
    L.append(f"Chains: {len(b)}; within 10% of best: {w} ({100*w//len(b)}%).")
except Exception:
    pass
os.makedirs("figs_and_results", exist_ok=True)
open("figs_and_results/RESULTS.md", "w").write("\n".join(L) + "\n")
print("wrote figs_and_results/RESULTS.md")
