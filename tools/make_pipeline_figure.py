#!/usr/bin/env python3
"""README figure: the automated SubMIT pipeline as four stage columns.
All connectors are horizontal/vertical segments joined by right-angle bends.
usage (repository root): python tools/make_pipeline_figure.py -> docs/figures/pipeline.png
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs", "figures")
INK, GREY = "#2b2b2b", "#6b6b6b"
STAGES = [  # title, accent, light fill
    ("1  Data", "#3b6fb6", "#eaf0f8"),
    ("2  Preparation", "#4d9a8a", "#e8f3f1"),
    ("3  Bayesian inversion", "#c0632b", "#fbeee6"),
    ("4  Results", "#7a5aa6", "#f0ebf6"),
]

fig = plt.figure(figsize=(15, 5.8))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 150); ax.set_ylim(16, 74); ax.axis("off")

# banner
ax.add_patch(FancyBboxPatch((3, 66.2), 144, 6, boxstyle="round,pad=0,rounding_size=1.2",
                            fc="#2b2b2b", ec="none"))
ax.text(75, 69.2, "Fully automated:  one command from downloaded waveforms to final figures  —  no manual tuning",
        ha="center", va="center", color="white", fontsize=14, fontweight="bold")

W, H = 30.0, 7.6                       # box size
COLX = [5, 42, 79, 116]                # left edge of each column's boxes
GAP = 3.5                              # vertical gap between boxes


def panel(i, top, bottom):
    t, acc, fill = STAGES[i]
    x = COLX[i] - 2.5
    ax.add_patch(FancyBboxPatch((x, bottom), W + 5, top - bottom, boxstyle="round,pad=0,rounding_size=1.5",
                                fc=fill, ec="none", zorder=0))
    ax.text(x + 1.5, top - 2.3, t, ha="left", va="center", fontsize=13, fontweight="bold", color=acc)


def box(i, y, text, bold=None):
    """box in column i with its top at y; returns (x0, y0, x1, y1)"""
    acc = STAGES[i][1]
    x0 = COLX[i]
    ax.add_patch(FancyBboxPatch((x0, y - H), W, H, boxstyle="round,pad=0,rounding_size=0.8",
                                fc="white", ec=acc, lw=1.6, zorder=2))
    if bold:
        ax.text(x0 + W / 2, y - 2.3, bold, ha="center", va="center", fontsize=10.5,
                fontweight="bold", color=INK, zorder=3)
        ax.text(x0 + W / 2, y - 5.2, text, ha="center", va="center", fontsize=9, color=GREY,
                linespacing=1.25, zorder=3)
    else:
        ax.text(x0 + W / 2, y - H / 2, text, ha="center", va="center", fontsize=10, color=INK, zorder=3)
    return (x0, y - H, x0 + W, y)


def path(pts, color=INK, lw=1.6, head=True):
    """polyline through right-angle corners, arrowhead on the last segment"""
    xs, ys = zip(*pts)
    ax.plot(xs[:-1] + (xs[-1],), ys[:-1] + (ys[-1],), color=color, lw=lw, zorder=1,
            solid_joinstyle="miter", solid_capstyle="butt")
    if head:
        (x0, y0), (x1, y1) = pts[-2], pts[-1]
        ax.annotate("", (x1, y1), (x0 + 0.8 * (x1 - x0), y0 + 0.8 * (y1 - y0)),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=16), zorder=4)


def down(b_from, b_to, color=INK):
    xm = (b_from[0] + b_from[2]) / 2
    path([(xm, b_from[1]), (xm, b_to[3])], color)


TOP = 60
# column 1 - data
c1 = [box(0, TOP, "teleseismic (30°–90°) + regional", "Download waveforms"),
      box(0, TOP - (H + GAP), "response removal, rotation,\narrival picks, 1 sample/s", "Step 1: processing")]
# column 2 - preparation
c2 = [box(1, TOP, "amplitude outliers, noise,\ndead channels, drift", "Screening, stage 1"),
      box(1, TOP - (H + GAP), "mtel3 (teleseismic), fk (regional),\nCRUST1.0 structure", "Green's functions"),
      box(1, TOP - 2 * (H + GAP), "aftershock map, catalog moment\ntensor, bounds, windows", "Priors & configuration")]
# column 3 - inversion
y = TOP
c3 = []
for txt, b in (("time bounds, windows, weights,\nscreening stage 2 (polarity, gain)", "N = 1, then calibrate"),
               ("N = 1 re-run with the calibrated\nsettings; N = 2, 3, ... (32 chains)", "MCMC inversions")):
    c3.append(box(2, y, txt, b)); y -= H + GAP
# decision diamond
dcx, dcy, dw, dh = COLX[2] + W / 2, y - 5.2, 15.5, 5.2
ax.add_patch(Polygon([(dcx, dcy + dh), (dcx + dw, dcy), (dcx, dcy - dh), (dcx - dw, dcy)],
                     closed=True, fc="white", ec=STAGES[2][1], lw=1.6, zorder=2))
ax.text(dcx, dcy + 0.9, "smallest N within 5% of", ha="center", va="center", fontsize=9.5, color=INK, zorder=3)
ax.text(dcx, dcy - 1.5, "the best fit, below max N?", ha="center", va="center", fontsize=9.5, color=INK, zorder=3)
# column 4 - results
c4 = [box(3, TOP, "selected model", "Forward modeling"),
      box(3, TOP - (H + GAP), "subevents, uncertainties, waveform\nfits, station maps, RESULTS.md", "Figures & summary")]

panel(0, TOP + 3.5, c1[-1][1] - 2.5)
panel(1, TOP + 3.5, c2[-1][1] - 2.5)
panel(2, TOP + 3.5, dcy - dh - 9.0)
panel(3, TOP + 3.5, c4[-1][1] - 2.5)
# redraw the boxes above the panels (panels were added afterwards)
for p in ax.patches:
    if isinstance(p, FancyBboxPatch) and p.get_facecolor()[:3] == (1, 1, 1):
        p.set_zorder(2)

# vertical arrows inside columns
down(c1[0], c1[1]); down(c2[0], c2[1]); down(c2[1], c2[2]); down(c3[0], c3[1]); down(c4[0], c4[1])
path([((c3[1][0] + c3[1][2]) / 2, c3[1][1]), (dcx, dcy + dh)])

# column-to-column connectors: out of the bottom box, through the gutter, into the next top box
def bridge(b_from, b_to):
    gx = (b_from[2] + b_to[0]) / 2
    ym_from = (b_from[1] + b_from[3]) / 2
    ym_to = (b_to[1] + b_to[3]) / 2
    path([(b_from[2], ym_from), (gx, ym_from), (gx, ym_to), (b_to[0], ym_to)])

bridge(c1[-1], c2[0])
bridge(c2[-1], c3[0])

# decision: yes -> results (right, up, right);  no -> add a subevent (down, left, up)
gx = (COLX[2] + W + COLX[3]) / 2
path([(dcx + dw, dcy), (gx, dcy), (gx, (c4[0][1] + c4[0][3]) / 2), (c4[0][0], (c4[0][1] + c4[0][3]) / 2)],
     color=STAGES[3][1])
ax.text(dcx + dw + 0.8, dcy + 1.2, "yes", fontsize=10, fontweight="bold", color=STAGES[3][1])
lx = COLX[2] - 0.8
ny = dcy - dh - 4.3
path([(dcx, dcy - dh), (dcx, ny), (lx, ny), (lx, (c3[1][1] + c3[1][3]) / 2),
      (c3[1][0], (c3[1][1] + c3[1][3]) / 2)], color=STAGES[2][1])
ax.text(dcx + 0.8, dcy - dh - 1.8, "no", fontsize=10, fontweight="bold", color=STAGES[2][1])
ax.text((dcx + lx) / 2 + 3, ny - 1.6, "add one subevent", ha="center", fontsize=9.5, color=STAGES[2][1])

fig.savefig(os.path.join(OUT, "pipeline.png"), dpi=150)
print("wrote docs/figures/pipeline.png")
