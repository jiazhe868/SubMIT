#!/usr/bin/env python3
"""Write figs_and_results/RESULTS.md: the selected model, the L-curve table (plot_lcurve.py), the
subevent table and the convergence report of the chains.

usage (in IRIS/, called by step4_figures.sh): python make_results_doc.py <event> <N>"""
import glob, os, sys
import numpy as np

ev, nn = sys.argv[1], int(sys.argv[2])
inv, fwd = f"inv_{ev}_{nn}sub", f"fwd_{ev}_{nn}sub"
model = np.loadtxt(f"{fwd}/Input.model", ndmin=2)
fm = np.loadtxt(f"{fwd}/fm.dat", ndmin=2)[:, 1:7]
def mw(m): return 2/3*np.log10(np.sqrt((m[0]**2+2*m[1]**2+2*m[2]**2+m[3]**2+2*m[4]**2+m[5]**2)/2)*1e27)-10.7
misfit = open(f"{fwd}/misfit.dat").read().strip()
# step 4's forward run re-evaluates the model with the current code. The two differ when the model was
# sampled with another code version: the published v2 examples in --mode figures (California 3sub:
# 1.7905 sampled, 1.8850 with the v3 timing of README section 4.2)
forward_misfit_note = ""
forward_log_path = f"{fwd}/ffwd.log"
if os.path.exists(forward_log_path):
    for log_line in open(forward_log_path):
        if log_line.strip().startswith("penalized misfit"):
            forward_misfit = float(log_line.split()[-1])
            if abs(forward_misfit / float(misfit) - 1) > 1e-3:
                forward_misfit_note = (f" (this code re-evaluates the model as {forward_misfit:.4f}: it was sampled "
                                       f"with another code version)")
mwcat = ""
try:
    mwcat = f" (catalog Mw {float(open(f'{inv}/mainshock.dat').read().split()[3]):.1f})"
except Exception:
    pass
L = [f"# Subevent inversion results: {ev}", "",
     f"- Selected subevents: **{nn}** | best misfit **{misfit}**{forward_misfit_note}",
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
      "- `lcurve.pdf` — model selection: best misfit per number of subevents (star: the selected number, the "
      "smallest of those tied with the lowest AIC; dashed line: a misfit 5% lower per added subevent, the "
      "cost of a subevent; hollow points: numbers that cannot be selected, their last subevent lies on the "
      "bound of the centroid range; small open circles: island bests) and, to its right, the factors of the misfit", ""]
if os.path.exists(f"{inv}/convergence.txt"):
    L += ["## Convergence of the chains", "```", open(f"{inv}/convergence.txt").read().rstrip(), "```", ""]
os.makedirs("figs_and_results", exist_ok=True)
open("figs_and_results/RESULTS.md", "w").write("\n".join(L) + "\n")
print("wrote figs_and_results/RESULTS.md")
