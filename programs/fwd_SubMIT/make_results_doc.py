#!/usr/bin/env python3
"""Write figs_and_results/RESULTS.md: the selected model, the L-curve table (plot_lcurve.py), the
subevent table and the convergence report of the chains.

The subevent table also lists what the rules of the search act on (finv.f90, forward_model.c), so a
reader does not have to work them out from the figure (added 2026-11-17, PITFALLS B4 and D6):
- the distance of every subevent from E1, which sits at the hypocentre, and that distance divided by
  the time since E1 began (its centroid time minus half its duration). The causality rule caps this
  speed. The speed between two centroids read off the figure is not capped and can be higher when
  the subevents overlap in time;
- the nodal planes of every subevent and the angle between its moment tensor and the double couple
  of the total tensor, the quantity that the mechanism prior acts on (a strike rotated by 10 degrees
  on a vertical strike-slip fault is an angle of 20 degrees).

usage (in IRIS/, called by step4_figures.sh): python make_results_doc.py <event> <N>"""
import glob, os, sys
import numpy as np

ev, nn = sys.argv[1], int(sys.argv[2])
inv, fwd = f"inv_{ev}_{nn}sub", f"fwd_{ev}_{nn}sub"
model = np.loadtxt(f"{fwd}/Input.model", ndmin=2)
fm = np.loadtxt(f"{fwd}/fm.dat", ndmin=2)[:, 1:7]
def mw(m): return 2/3*np.log10(np.sqrt((m[0]**2+2*m[1]**2+2*m[2]**2+m[3]**2+2*m[4]**2+m[5]**2)/2)*1e27)-10.7


def tensor_matrix(six_components):
    """3 x 3 moment tensor (x north, y east, z down) from the fm.dat order Mxx Mxy Mxz Myy Myz Mzz"""
    mxx, mxy, mxz, myy, myz, mzz = six_components
    return np.array([[mxx, mxy, mxz], [mxy, myy, myz], [mxz, myz, mzz]])


def double_couple_axes(matrix):
    """unit pressure and tension axes of a moment tensor (eigenvectors of its smallest and largest eigenvalue)"""
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return eigenvectors[:, 0], eigenvectors[:, 2]


def nodal_planes(matrix):
    """(strike, dip, rake) in degrees of the two nodal planes of the best double couple (Aki & Richards convention)"""
    pressure_axis, tension_axis = double_couple_axes(matrix)
    planes = []
    for sign in (1, -1):
        normal = (tension_axis + sign * pressure_axis) / np.sqrt(2)
        slip = (tension_axis - sign * pressure_axis) / np.sqrt(2)
        if normal[2] > 0:                       # the normal of the hanging wall points up (z is down)
            normal, slip = -normal, -slip
        strike = np.degrees(np.arctan2(-normal[0], normal[1])) % 360
        dip = np.degrees(np.arccos(-normal[2]))
        along_strike = np.array([np.cos(np.radians(strike)), np.sin(np.radians(strike)), 0.0])
        rake = np.degrees(np.arctan2(-slip[2] / max(np.sin(np.radians(dip)), 1e-9), slip @ along_strike))
        planes.append((strike, dip, rake))
    return planes


def angle_between_tensors(first, second):
    """angle in degrees between two moment tensors as vectors of their nine components (0-180)"""
    cosine = np.sum(first * second) / np.sqrt(np.sum(first * first) * np.sum(second * second))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))
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
# the reference of the mechanism prior: the double couple of the catalogue tensor (totalmt.dat: Mxx Myy Mxy Mxz Myz)
reference_double_couple = None
try:
    mxx, myy, mxy, mxz, myz = np.loadtxt(f"{inv}/totalmt.dat").ravel()[:5]
    pressure_axis, tension_axis = double_couple_axes(tensor_matrix((mxx, mxy, mxz, myy, myz, -mxx - myy)))
    reference_double_couple = np.outer(tension_axis, tension_axis) - np.outer(pressure_axis, pressure_axis)
except Exception:
    pass
first_begins = model[0][0] - model[0][3] / 2          # the time at which the rupture front leaves E1
L += ["## Best model", "",
      "| # | centroid (s) | X east (km) | Y north (km) | duration (s) | depth (km) | Mw | from E1 (km) | speed since E1 began (km/s) "
      "| nodal planes (strike/dip/rake) | angle from the total mechanism (deg) |",
      "|---|---|---|---|---|---|---|---|---|---|---|"]
for i, r in enumerate(model):
    distance_from_first = float(np.sqrt((r[1] - model[0][1])**2 + (r[2] - model[0][2])**2 + (r[6] - model[0][6])**2))
    speed = f"{distance_from_first / (r[0] - first_begins):.2f}" if i > 0 and r[0] > first_begins else "-"
    matrix = tensor_matrix(fm[i])
    planes = " ; ".join(f"{strike:.0f}/{dip:.0f}/{rake:.0f}" for strike, dip, rake in nodal_planes(matrix))
    angle = f"{angle_between_tensors(matrix, reference_double_couple):.0f}" if reference_double_couple is not None else "-"
    L.append(f"| E{i+1} | {r[0]:.2f} | {r[1]:.2f} | {r[2]:.2f} | {r[3]:.2f} | {r[6]:.1f} | {mw(fm[i]):.2f} "
             f"| {distance_from_first:.0f} | {speed} | {planes} | {angle} |")
L += ["", "E1 sits at the hypocentre when there are two or more subevents. The speed is the distance from E1 divided by "
          "the time since E1 began (its centroid time minus half its duration); the search caps it (`vmax_cau` in "
          "`search_par.file`). The angle is measured between each moment tensor and the double couple of the catalogue "
          "tensor; the mechanism prior acts on it (`MechanismSigma` in `Par.file`)."]
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
