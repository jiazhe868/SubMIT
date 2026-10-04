#!/usr/bin/env python3
"""GMT-free waveform-fit plates for SubMIT step 4 (matplotlib + obspy only).

Run inside a fwd_* directory after ./ffwd. Reads Par.file, station info files
and waveforms/{P,SH,rayl?}_{obs,syn}_%04d.sac; writes one multi-page PDF per
wave type: fits_P.pdf, fits_Pvel.pdf, fits_SH.pdf, fits_rayl.pdf.
Replaces plotP/plotPvel/plotSH/plotrayl GMT4+pssac2 scripts.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from obspy import read

PER_COL, NCOL = 12, 3   # 36/page: tighter plates, fewer pages (user 2026-08-13)  # stations per column, columns per page
PLOT_DT = 0.2          # traces are interpolated to this dt for smooth plotting
LW_OBS, LW_SYN = 2.0, 2.0

def smooth(tr):
    """weighted-average-slopes interpolation for display smoothness"""
    try:
        tr.interpolate(sampling_rate=1.0/PLOT_DT, method="weighted_average_slopes")
    except Exception:
        pass
    return tr

def par(key, path="Par.file"):
    for line in open(path):
        if line.split("=")[0].strip() == key.rstrip("="):
            return float(line.split("=")[1].split("#")[0])
    raise KeyError(key)

def stations(path, three_comp=False):
    out = []
    for line in open(path):
        t = line.split()
        if not t:
            continue
        name = ".".join(os.path.basename(t[0]).split(".")[:2])  # NET.STA only
        if three_comp:
            out.append((name, float(t[1]), float(t[2])))   # name, stlo, stla
        else:
            out.append((name, float(t[1]), float(t[2])))   # name, dist(deg), az
    return out

def plate(pdfname, pairs, bg, nd, title, cut=0.0):
    """pairs: list of (label, sublabel, obsfile, synfile) already sorted.
    cut: seconds hidden at each end (the maxshft alignment slack for this
    wave type) - traces are drawn in full but the window is truncated."""
    if not pairs:
        return
    # BALANCED pagination: never a nearly-empty trailing page (56 traces ->
    # 19/19/18, not 27/27/2); per-column count shrinks accordingly
    per_page = PER_COL * NCOL
    npages = max(1, -(-len(pairs) // per_page))
    per_page = -(-len(pairs) // npages)
    pages = [pairs[i:i + per_page] for i in range(0, len(pairs), per_page)]
    with PdfPages(pdfname) as pdf:
        for ip, chunk in enumerate(pages):
            fig, axes = plt.subplots(1, NCOL, figsize=(8.5, 11), sharey=False)
            axes = np.atleast_1d(axes)
            for ic in range(NCOL):
                ax = axes[ic]
                percol = -(-per_page // NCOL)
                sub = chunk[ic * percol:(ic + 1) * percol]
                for k, (lab, sub2, fo, fs) in enumerate(sub):
                    y0 = len(sub) - k
                    try:
                        o = smooth(read(fo)[0]); s = smooth(read(fs)[0])
                    except Exception:
                        continue
                    t = bg + np.arange(o.stats.npts) * o.stats.delta
                    ts = bg + np.arange(s.stats.npts) * s.stats.delta
                    # normalize on the VISIBLE window only
                    vo = o.data[(t >= bg + cut) & (t <= nd - cut)]
                    vs = s.data[(ts >= bg + cut) & (ts <= nd - cut)]
                    a = max(np.max(np.abs(vo), initial=0), np.max(np.abs(vs), initial=0), 1e-30)
                    ax.plot(t, o.data / a * 0.45 + y0, color="0.1", lw=LW_OBS)
                    ax.plot(ts, s.data / a * 0.45 + y0, color="crimson", lw=LW_SYN)
                    ax.text(bg + cut, y0 + 0.26, lab, fontsize=10.5, va="bottom", fontweight="bold")
                    ax.text(bg + cut, y0 - 0.44, sub2, fontsize=9.5, va="bottom", color="0.25")
                ax.set_xlim(bg + cut, nd - cut)
                ax.set_ylim(0.2, -(-per_page // NCOL) + 1)
                ax.set_yticks([])
                ax.set_xlabel("time (s)", fontsize=10)
                for sp in ("top", "right", "left"):
                    ax.spines[sp].set_visible(False)
                ax.tick_params(labelsize=9)
            fig.suptitle(f"{title} \N{EM DASH} obs (black) vs syn (red)"
                         + (f"  ({ip+1}/{len(pages)})" if len(pages) > 1 else ""),
                         fontsize=11)
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig)
            plt.close(fig)
    print(f"wrote {pdfname} ({len(pairs)} traces, {len(pages)} page(s))")

def par_or(key, default):
    try:
        return par(key)
    except KeyError:
        return default

def main():
    w = "waveforms"
    cutP = par_or("maxshftP", 2.0)
    cutSH = par_or("maxshftSH", 5.0)
    cutR = par_or("maxshftRayl", 3.0)
    # --- teleseismic P: displacement then vel_ (same index sequence)
    staP = stations("stations.info")
    ndisp = sum(1 for n, _, _ in staP if not n.startswith("vel_"))
    bg, nd = par("bg_timeP"), par("nd_timeP")
    def ppair(i, n, d, a):
        return (n.replace("vel_", ""), f"{d:.1f}\N{DEGREE SIGN} az {a:.0f}\N{DEGREE SIGN}",
                f"{w}/P_obs_{i:04d}.sac", f"{w}/P_syn_{i:04d}.sac")
    disp = [(a, ppair(i, n, d, a)) for i, (n, d, a) in enumerate(staP) if i < ndisp]
    vel = [(a, ppair(i, n, d, a)) for i, (n, d, a) in enumerate(staP) if i >= ndisp]
    plate("fits_P.pdf", [p for _, p in sorted(disp)], bg, nd, "Teleseismic P (displacement)", cut=cutP)
    plate("fits_Pvel.pdf", [p for _, p in sorted(vel)], bg, nd, "Teleseismic P (velocity)", cut=cutP)
    # --- SH
    staSH = stations("stationsSH.info")
    bg, nd = par("bg_timeSH"), par("nd_timeSH")
    sh = sorted((a, (n, f"{d:.1f}\N{DEGREE SIGN} az {a:.0f}\N{DEGREE SIGN}",
                     f"{w}/SH_obs_{i:04d}.sac", f"{w}/SH_syn_{i:04d}.sac"))
                for i, (n, d, a) in enumerate(staSH))
    plate("fits_SH.pdf", [p for _, p in sh], bg, nd, "Teleseismic SH", cut=cutSH)
    # --- regional 3-component: one ROW per station, E/N/Z as the three columns,
    # labeled with epicentral distance (km) instead of lon/lat
    staL = stations("stationsloc.info", three_comp=True)
    evlo, evla = par("evlo"), par("evla")
    bg, nd = par("bg_timeRayl"), par("nd_timeRayl")
    def dist_km(lo, la):
        return np.hypot((lo - evlo) * 111.32 * np.cos(np.radians(evla)),
                        (la - evla) * 110.574)
    rows = sorted((dist_km(lo, la), i, n) for i, (n, lo, la) in enumerate(staL))
    per_page = PER_COL
    pages = [rows[i:i + per_page] for i in range(0, len(rows), per_page)]
    with PdfPages("fits_rayl.pdf") as pdf:
        for ip, chunk in enumerate(pages):
            fig, axes = plt.subplots(1, 3, figsize=(8.5, 11))
            for ic, comp in enumerate("enz"):
                ax = axes[ic]
                for k, (dk, i, n) in enumerate(chunk):
                    y0 = len(chunk) - k
                    try:
                        o = smooth(read(f"{w}/rayl{comp}_obs_{i:04d}.sac")[0])
                        sy = smooth(read(f"{w}/rayl{comp}_syn_{i:04d}.sac")[0])
                    except Exception:
                        continue
                    t = bg + np.arange(o.stats.npts) * o.stats.delta
                    ts = bg + np.arange(sy.stats.npts) * sy.stats.delta
                    vo = o.data[(t >= bg + cutR) & (t <= nd - cutR)]
                    vs = sy.data[(ts >= bg + cutR) & (ts <= nd - cutR)]
                    a = max(np.max(np.abs(vo), initial=0), np.max(np.abs(vs), initial=0), 1e-30)
                    ax.plot(t, o.data / a * 0.45 + y0, color="0.1", lw=LW_OBS)
                    ax.plot(ts, sy.data / a * 0.45 + y0, color="crimson", lw=LW_SYN)
                    if ic == 0:
                        ax.text(bg + cutR, y0 + 0.28, f"{n}  {dk:.0f} km", fontsize=10.5, va="bottom", fontweight="bold")
                ax.set_title(comp.upper(), fontsize=10)
                ax.set_xlim(bg + cutR, nd - cutR); ax.set_ylim(0.2, per_page + 1)
                ax.set_yticks([]); ax.set_xlabel("time (s)", fontsize=10)
                for sp in ("top", "right", "left"):
                    ax.spines[sp].set_visible(False)
                ax.tick_params(labelsize=9)
            fig.suptitle("Regional full waveforms — obs (black) vs syn (red)"
                         + (f"  ({ip+1}/{len(pages)})" if len(pages) > 1 else ""), fontsize=11)
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig); plt.close(fig)
    print(f"wrote fits_rayl.pdf ({len(rows)} stations x 3 comps, {len(pages)} page(s))")

if __name__ == "__main__":
    main()
