#!/usr/bin/env python3
"""GMT-free subevent model figure for SubMIT step 4 (matplotlib + obspy only).

Run inside a fwd_* directory after ./ffwd. Reads Input.model (cen x y dura vr
theta depth per row), fm.dat (1e27 Mxx Mxy Mxz Myy Myz Mzz per subevent) and
stf/stf_%04d.sac + stf/mrf_sum.sac. Writes subevents.pdf with three panels:
map view with beachballs, depth cross-section, and source time functions.
Replaces subeveplot/plot.sh (GMT4 psmeca/pssac2).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from obspy import read
from obspy.imaging.beachball import beach

def mw_of(m6):
    m0 = np.sqrt((m6[0]**2 + 2*m6[1]**2 + 2*m6[2]**2 + m6[3]**2 + 2*m6[4]**2 + m6[5]**2) / 2) * 1e27
    return 2.0/3.0 * np.log10(m0) - 10.7

def xyz2rtp(m6):
    """(Mxx,Mxy,Mxz,Myy,Myz,Mzz), x north y east z down -> (rr,tt,pp,rt,rp,tp)"""
    mxx, mxy, mxz, myy, myz, mzz = m6
    return [mzz, mxx, myy, mxz, -myz, -mxy]

def section_rtp(m6):
    """Vertical-plane beachball, viewed from the SOUTH looking north.
    Fake geographic frame: fake-down = +north (projection axis), fake-north =
    up (-z, so 'up' in the plot is shallower), fake-east = east. Rotating
    M (x=n, y=e, z=d) into that frame and converting to USE gives:
    (rr,tt,pp,rt,rp,tp) = (Mxx, Mzz, Myy, -Mxz, -Mxy, +Myz)."""
    mxx, mxy, mxz, myy, myz, mzz = m6
    # E-W mirror applied (calibrated against a known east-dipping plane: the
    # unmirrored rendering appeared west-dipping): flips the two odd-in-east
    # components relative to the plain rotation
    return [mxx, mzz, myy, -mxz, mxy, -myz]

model = np.loadtxt("Input.model", ndmin=2)   # cen x y dura vr theta depth
fm = np.loadtxt("fm.dat", ndmin=2)[:, 1:7]
nsub = len(model)
mws = [mw_of(fm[i]) for i in range(nsub)]
tot = fm.sum(axis=0)
colors = plt.cm.tab10(np.linspace(0, 1, 10))

plt.rcParams.update({"font.size": 12})
fig = plt.figure(figsize=(14, 8))
# map on the left; moment-rate functions, total mechanism and depth section
# in separate panels on the right (no insets overlapping the map axes)
axm = fig.add_axes([0.07, 0.08, 0.50, 0.82])   # map view
axs = fig.add_axes([0.64, 0.60, 0.21, 0.30])   # moment-rate functions
axt = fig.add_axes([0.87, 0.60, 0.11, 0.30])   # total moment tensor
axz = fig.add_axes([0.64, 0.08, 0.34, 0.42])   # depth section

x, y, z = model[:, 1], model[:, 2], model[:, 6]
# hypocenter for lon/lat axes
evlo = evla = None
try:
    for line in open("Par.file"):
        k = line.split("=")[0].strip()
        if k == "evlo": evlo = float(line.split("=")[1].split("#")[0])
        if k == "evla": evla = float(line.split("=")[1].split("#")[0])
except Exception:
    pass
def km2lon(xk): return evlo + xk / (111.32 * np.cos(np.radians(evla)))
def km2lat(yk): return evla + yk / 110.574
# 95% location intervals (2.5-97.5 percentiles) from the MCMC ensemble
perc = None
try:
    smp = np.loadtxt("allsamples.dat", ndmin=2)
    perc = {}
    for i in range(nsub):
        for nmch, col in (("x", 2 + 5*i), ("y", 3 + 5*i), ("z", 5 + 5*i)):
            lo, hi = np.percentile(smp[:, col], [2.5, 97.5])
            perc[(i, nmch)] = (lo, hi)
except Exception:
    perc = None
pad = max(20.0, 0.4 * max(x.ptp(), y.ptp(), 1))
axm.set_xlim(km2lon(min(x.min(), 0) - pad), km2lon(max(x.max(), 0) + pad))
axm.set_ylim(km2lat(min(y.min(), 0) - pad), km2lat(max(y.max(), 0) + pad))
axm.set_aspect(1.0 / np.cos(np.radians(evla)))   # isotropic km on a deg grid
axm.set_anchor("E")   # hug the right-hand panels when the aspect narrows the map
sizes = [0.12 * (mw / max(mws)) ** 3 * pad for mw in mws]
deg_x = lambda km: km / (111.32 * np.cos(np.radians(evla)))
deg_y = lambda km: km / 110.574
for i in range(nsub):
    if perc:
        (xl, xh), (yl, yh) = perc[(i, "x")], perc[(i, "y")]
        axm.errorbar(km2lon(x[i]), km2lat(y[i]),
                     xerr=[[deg_x(max(x[i]-xl, 0))], [deg_x(max(xh-x[i], 0))]],
                     yerr=[[deg_y(max(y[i]-yl, 0))], [deg_y(max(yh-y[i], 0))]],
                     fmt="none", ecolor=colors[i % 10], elinewidth=0.9,
                     capsize=2.5, zorder=7)
    # width as (x,y) tuple: with the 1/cos(lat) aspect a scalar (lon-deg)
    # width renders as an ellipse; matched deg widths give a screen circle
    bw = (2*deg_x(sizes[i]), 2*deg_y(sizes[i]))
    b = beach(xyz2rtp(fm[i]), xy=(km2lon(x[i]), km2lat(y[i])),
              width=bw, facecolor=colors[i % 10], linewidth=0.6)
    b.set_zorder(5)
    axm.add_collection(b)
    # best-double-couple nodal planes on top of the full-MT fill
    try:
        from obspy.imaging.beachball import MomentTensor, mt2plane
        npl = mt2plane(MomentTensor(xyz2rtp(fm[i]), 0))
        bn = beach((npl.strike, npl.dip, npl.rake), xy=(km2lon(x[i]), km2lat(y[i])),
                   width=bw, nofill=True, linewidth=0.8, edgecolor="k")
        bn.set_zorder(6)
        axm.add_collection(bn)
    except Exception:
        pass
axm.plot(evlo, evla, "k*", ms=16, zorder=6)
# labels: try positions around each beachball (above, below, right, left,
# diagonals) and keep the first whose box overlaps no label placed so far
fig.canvas.draw()
from matplotlib.transforms import Bbox
balls = []   # labels must not cover any beachball (other than touching their own)
for j in range(nsub):
    (x0, y0), (x1, y1) = axm.transData.transform(
        [(km2lon(x[j] - sizes[j]), km2lat(y[j] - sizes[j])),
         (km2lon(x[j] + sizes[j]), km2lat(y[j] + sizes[j]))])
    balls.append(Bbox([[x0, y0], [x1, y1]]))
placed = []
for i in sorted(range(nsub), key=lambda k: -mws[k]):
    txt = f"E{i+1} Mw{mws[i]:.1f}  {model[i,0]:.1f}s"
    r = sizes[i] + 0.04 * pad
    cands = [(0, r, "center", "bottom"), (0, -r, "center", "top"),
             (r, 0, "left", "center"), (-r, 0, "right", "center"),
             (r, r, "left", "bottom"), (-r, r, "right", "bottom"),
             (r, -r, "left", "top"), (-r, -r, "right", "top")]
    for k, (dx, dy, ha, va) in enumerate(cands + [(0, 2.6 * r, "center", "bottom")]):
        a = axm.annotate(txt, (km2lon(x[i] + dx), km2lat(y[i] + dy)), ha=ha, va=va,
                         fontsize=11, fontweight="bold", zorder=10,
                         bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.0))
        bb = a.get_window_extent(fig.canvas.get_renderer()).expanded(1.05, 1.15)
        hit = any(bb.overlaps(q) for q in placed) or \
            any(bb.overlaps(balls[j]) for j in range(nsub) if j != i)
        if not hit or k == len(cands):
            placed.append(bb)
            break
        a.remove()

axt.set_xlim(0, 1); axt.set_ylim(0, 1); axt.set_aspect("equal"); axt.axis("off")
bt = beach(xyz2rtp(tot), xy=(0.5, 0.5), width=0.8, facecolor="0.4", linewidth=0.6)
axt.add_collection(bt)
try:
    from obspy.imaging.beachball import MomentTensor, mt2plane
    npl = mt2plane(MomentTensor(xyz2rtp(tot), 0))
    axt.add_collection(beach((npl.strike, npl.dip, npl.rake), xy=(0.5, 0.5), width=0.8,
                             nofill=True, linewidth=0.8, edgecolor="k"))
except Exception:
    pass
axt.set_title(f"Total Mw {mw_of(tot):.2f}", fontsize=11)
# distance scale bar (nice round length ~1/4 of the x span)
sspan = (max(x.max(), 0) - min(x.min(), 0) + 2 * pad)
sbar = min([1, 2, 5, 10, 20, 50, 100, 200], key=lambda v: abs(v - sspan / 4.0))
sx1 = km2lon(max(x.max(), 0) + pad - sbar - 0.06 * sspan)
sy = km2lat(min(y.min(), 0) - pad + 0.07 * sspan)
axm.plot([sx1, sx1 + deg_x(sbar)], [sy, sy], "k-", lw=2.5, zorder=10,
         solid_capstyle="butt")
axm.annotate(f"{sbar} km", (sx1 + deg_x(sbar / 2.0), sy), ha="center",
             va="bottom", fontsize=10, zorder=10, xytext=(0, 3),
             textcoords="offset points")
axm.set_xlabel("Longitude (°)"); axm.set_ylabel("Latitude (°)")
axm.grid(alpha=0.25)

# depth cross-section (X-Z): side-view beachball projection, not the map view
xspan_s = max(x.max(), 0) - min(x.min(), 0) + 2 * pad
zpad = max(8.0, 0.5 * max(z.ptp(), 1))
zspan_s = z.ptp() + 2 * zpad
ve = max(1.0, round(xspan_s / zspan_s / 1.6))   # vertical exaggeration
ssec = [min(sz, 0.22 * zspan_s * ve) for sz in sizes]   # x-radius (km)
for i in range(nsub):
    if perc:
        (xl, xh), (zl, zh) = perc[(i, "x")], perc[(i, "z")]
        axz.errorbar(x[i], z[i],
                     xerr=[[max(x[i]-xl, 0)], [max(xh-x[i], 0)]],
                     yerr=[[max(z[i]-zl, 0)], [max(zh-z[i], 0)]],
                     fmt="none", ecolor=colors[i % 10], elinewidth=0.9,
                     capsize=2.5, zorder=7)
    b = beach(section_rtp(fm[i]), xy=(x[i], z[i]), width=(2*ssec[i], 2*ssec[i]/ve),
              facecolor=colors[i % 10], linewidth=0.6)
    b.set_zorder(5)
    axz.add_collection(b)
    try:
        from obspy.imaging.beachball import MomentTensor, mt2plane
        npl = mt2plane(MomentTensor(section_rtp(fm[i]), 0))
        bn = beach((npl.strike, npl.dip, npl.rake), xy=(x[i], z[i]),
                   width=(2*ssec[i], 2*ssec[i]/ve), nofill=True, linewidth=0.8, edgecolor="k")
        bn.set_zorder(6)
        axz.add_collection(bn)
    except Exception:
        pass
    axz.annotate(f"E{i+1}", (x[i], z[i] - ssec[i] / ve - 0.5), ha="center", va="bottom",
                 fontsize=11, fontweight="bold")
axz.set_xlim(min(x.min(), 0) - pad, max(x.max(), 0) + pad)   # km (map is deg)
axz.set_ylim(z.max() + zpad, max(z.min() - zpad, -2))   # depth increases downward
axz.set_aspect(ve)
axz.set_xlabel("X, east of hypocenter (km)"); axz.set_ylabel("depth (km)")
axz.set_title("Depth section (view from the south)"
              + (f", vertical exaggeration ×{ve:.0f}" if ve > 1 else ""), fontsize=11)
axz.grid(alpha=0.25)

# source time functions: stf_*.sac are unit-shape MRFs - scale each by its
# subevent scalar moment (from fm.dat) so relative amplitudes are physical;
# time axis shifted by stf_btime so the rupture starts near t=0
m0s = [np.sqrt((fm[i][0]**2 + 2*fm[i][1]**2 + 2*fm[i][2]**2 + fm[i][3]**2
                + 2*fm[i][4]**2 + fm[i][5]**2) / 2) for i in range(nsub)]
stf_bt = 0.0
try:
    for line in open("Par.file"):
        if line.split("=")[0].strip() == "stf_btime":
            stf_bt = float(line.split("=")[1].split("#")[0])
except Exception:
    pass
tsum, dsum = None, None
for i in range(nsub):
    tr = read(f"stf/stf_{i:04d}.sac")[0]
    d = tr.data / max(np.trapz(np.abs(tr.data)), 1e-30) * m0s[i]
    t = stf_bt + np.arange(tr.stats.npts) * tr.stats.delta
    axs.fill_between(t, 0, d, color=colors[i % 10], alpha=0.55, lw=0,
                     label=f"E{i+1}")
    if dsum is None:
        tsum, dsum = t, d.copy()
    else:
        n = min(len(dsum), len(d)); dsum = dsum[:n] + d[:n]; tsum = tsum[:n]
if dsum is not None:
    axs.plot(tsum, dsum, "k", lw=1.2, label="total")
axs.legend(fontsize=9, frameon=False, ncol=1, loc="upper right")
if dsum is not None:   # trim the flat tail so the pulses fill the panel
    act = np.where(dsum > 0.01 * dsum.max())[0]
    t_end = tsum[act[-1]] if len(act) else tsum[-1]
    axs.set_xlim(max(stf_bt, -2), t_end + 0.15 * (t_end - max(stf_bt, -2)) + 2)
axs.set_ylim(0, None)
axs.set_xlabel("time after origin (s)", fontsize=10)
axs.set_ylabel("moment rate", fontsize=10)
axs.set_yticks([]); axs.tick_params(labelsize=9)
axs.set_title("Moment-rate functions", fontsize=11)

fig.suptitle(f"Subevent model ({nsub} subevents), penalized misfit "
             f"{open('misfit.dat').read().strip() if __import__('os').path.exists('misfit.dat') else '?'}",
             fontsize=12)
with PdfPages("subevents.pdf") as pdf:
    pdf.savefig(fig)
fig.savefig("subevents.png", dpi=140)
print(f"wrote subevents.pdf / .png ({nsub} subevents, total Mw{mw_of(tot):.2f})")
import os
if os.path.exists("mainshock.dat"):
    try:
        mwcat = float(open("mainshock.dat").read().split()[3])
        if abs(mw_of(tot) - mwcat) > 0.2:
            print(f"WARNING: total Mw {mw_of(tot):.2f} differs from catalog Mw {mwcat:.1f} "
                  f"by >0.2 - check Tikhonov/DCconstrain/depth range")
    except Exception:
        pass
