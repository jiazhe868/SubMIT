#!/usr/bin/env python3
"""README figure: raw records (left) vs what the inversion uses (right).
Uses the screened data of --mode full runs of california and chile
(rejected records are moved to data/excluded/ by the screening). Run from the
repository root after those runs, or point SUBMIT_FIG_CA/_CH at the
inversion directories. Writes docs/figures/data_processing.png.
"""
import struct, glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

CA  = os.environ.get("SUBMIT_FIG_CA", "work/california/IRIS/inv_2024-12-05-mww70-off-coast-of-northern_3sub")
CH  = os.environ.get("SUBMIT_FIG_CH", "work/chile/IRIS/inv_2024-07-19-mww74-chile-argentina-border-region_4sub")

BLUE, RED, GRAY = "#3b6fb6", "#c0392b", "#9a9a9a"
C_S, C_R = "#e8853d", "#4d9a8a"

def sac(fn):
    b = open(fn, "rb").read()
    f = lambda w: struct.unpack("<f", b[4*w:4*w+4])[0]
    delta, bt, t1 = f(0), f(5), f(11)
    n = struct.unpack("<i", b[316:320])[0]
    d = np.array(struct.unpack(f"<{n}f", b[632:632+4*n]))
    return bt + delta*np.arange(n), d, t1

def norm(d): return d/np.abs(d).max()

plt.rcParams.update({"font.size": 11})
fig = plt.figure(figsize=(6.2, 9.2))
gs = GridSpec(5, 2, figure=fig, hspace=0.55, wspace=0.16,
              left=0.045, right=0.975, top=0.905, bottom=0.035)

def bare(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

def rowlabel(axL, txt):
    axL.text(0.0, 1.14, txt, transform=axL.transAxes, fontsize=11.5,
             fontweight="bold", color="#333333")

# column headers + flow arrow
fig.text(0.26, 0.955, "as recorded", ha="center", fontsize=14, color=RED, fontweight="bold")
fig.text(0.75, 0.955, "ready for the inversion", ha="center", fontsize=14, color=BLUE, fontweight="bold")
fig.text(0.505, 0.955, "→", ha="center", fontsize=17, color="#333333")

# ---- row 1: instrument fingerprint / drift removed --------------------------
axL = fig.add_subplot(gs[0, 0]); axR = fig.add_subplot(gs[0, 1])
t, d, _ = sac(f"{CA}/data/excluded/G.DZM.00.z")
dn = norm(d)
k = 81
trend = np.convolve(dn, np.ones(k)/k, mode="same")
axL.plot(t, dn, color=RED, lw=0.8)
clean = (dn - trend)[k:-k]
axR.plot(t[k:-k], norm(clean), color=BLUE, lw=0.8)
rowlabel(axL, "1  sensor fingerprint & drift removed")
bare(axL); bare(axR)

# ---- row 2: aligned on the wave arrival -------------------------------------
axL = fig.add_subplot(gs[1, 0]); axR = fig.add_subplot(gs[1, 1])
sta = []
with open(f"{CH}/stations.info") as f:
    for ln in f:
        p = ln.split()
        if p and p[0].endswith(".z"):
            sta.append((float(p[1]), p[0]))
sta.sort()
pick = [sta[3], sta[len(sta)//2], sta[-4]]
for i, (dist, name) in enumerate(pick):
    t, d, t1 = sac(f"{CH}/{name}")
    w = (t >= t1-15) & (t <= t1+80)
    m = np.abs(d[w]).max()
    axL.plot(t, d/m - 2.6*i, color=RED, lw=0.7, alpha=0.9)
    axR.plot(t[w]-t1, d[w]/m - 2.6*i, color=BLUE, lw=0.9)
axR.axvline(0, color="#333333", lw=0.9, ls="--")
rowlabel(axL, "2  aligned on the wave's arrival")
bare(axL); bare(axR)

# ---- row 3: impossible amplitudes rejected ----------------------------------
axL = fig.add_subplot(gs[2, 0]); axR = fig.add_subplot(gs[2, 1])
tp, dp, t1p = sac(f"{CH}/data/excluded/IM.PD31..z")
th, dh, t1h = sac(f"{CH}/data/BK.CMB.00.z")
tm, dm, t1m = sac(f"{CH}/data/excluded/AU.MAW..z")
ref = np.abs(dp).max()          # common absolute scale
for i, (tt, dd, tt1, c) in enumerate([(tp, dp, t1p, RED), (th, dh, t1h, BLUE), (tm, dm, t1m, RED)]):
    w = (tt >= tt1-20) & (tt <= tt1+160)
    axL.plot(tt[w]-tt1, dd[w]/ref - 2.4*i, color=c, lw=0.8)
def peak(tt, dd, tt1):
    w = (tt >= tt1-20) & (tt <= tt1+160)
    return np.abs(dd[w]).max()
def ratio_label(r):
    return f"×{r:.0f}" if r >= 1.5 else f"×{r:.2g}"
ph = peak(th, dh, t1h)   # amplitudes relative to the healthy station
axL.text(160, 0.9, ratio_label(peak(tp, dp, t1p) / ph), color=RED, fontsize=10, ha="right")
axL.text(160, -4.1, ratio_label(peak(tm, dm, t1m) / ph), color=RED, fontsize=10, ha="right")
w = (th >= t1h-20) & (th <= t1h+160)
axR.plot(th[w]-t1h, norm(dh[w]) - 2.4, color=BLUE, lw=0.9)
axL.set_ylim(-5.9, 1.6); axR.set_ylim(-5.9, 1.6)
rowlabel(axL, "3  impossible amplitudes rejected")
bare(axL); bare(axR)

# ---- row 4: broken channels rejected -----------------------------------------
axL = fig.add_subplot(gs[3, 0]); axR = fig.add_subplot(gs[3, 1])
t2, d2, t12 = sac(f"{CH}/data/excluded/GT.DBIC.00.t")
t3, d3, t13 = sac(f"{CH}/data/BK.CMB.00.t")
w2 = (t2 >= t12-20) & (t2 <= t12+180)
w3 = (t3 >= t13-20) & (t3 <= t13+180)
axL.plot(t2[w2]-t12, norm(d2[w2]) + 1.3, color=RED, lw=0.9)
axL.plot(t3[w3]-t13, norm(d3[w3]) - 1.3, color=BLUE, lw=0.9)
axL.text(175, 2.15, "upside-down", color=RED, fontsize=10, ha="right")
axR.plot(t3[w3]-t13, norm(d3[w3]) - 1.3, color=BLUE, lw=0.9)
axL.set_ylim(-2.9, 2.9); axR.set_ylim(-2.9, 2.9)
rowlabel(axL, "4  broken channels rejected")
bare(axL); bare(axR)

# ---- row 5: wave types weighted ----------------------------------------------
axL = fig.add_subplot(gs[4, 0]); axR = fig.add_subplot(gs[4, 1])
# Chile 4-subevent model: wave-type weights (Par.file) and weighted band
# residuals printed by ffwd ("bandresid:" line)
wgt = np.array([3.64624, 1.0, 0.050148]); rw = np.array([1.265707e-5, 9.484927e-6, 5.144819e-6])
raw = rw/wgt**2; raw = raw/raw.sum()*100
bal = rw/rw.sum()*100
cols = [BLUE, C_S, C_R]
names = ["P", "S", "surface"]
for ax, sh in [(axL, raw), (axR, bal)]:
    left = 0
    for s, c, nm in zip(sh, cols, names):
        ax.barh(0, s, left=left, color=c, height=0.62)
        if s > 20:
            ax.text(left+s/2, 0, f"{nm} {s:.0f}%", ha="center", va="center",
                    color="white", fontsize=10, fontweight="bold")
        elif s > 7:
            ax.text(left+s/2, 0, f"{s:.0f}%", ha="center", va="center",
                    color="white", fontsize=9, fontweight="bold")
        left += s
    ax.set_xlim(0, 100); ax.set_ylim(-1.1, 1.1)
    bare(ax)
axL.text(50, -0.92, "one wave type drowns the rest", ha="center", fontsize=9, color=RED)
axR.text(50, -0.92, "weighted by reliability", ha="center", fontsize=9, color=BLUE)
rowlabel(axL, "5  wave types weighted, not shouted")

out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs", "figures")
fig.savefig(os.path.join(out, "data_processing.png"), dpi=200)
print("wrote docs/figures/data_processing.png")
