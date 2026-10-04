#!/usr/bin/env python3
"""Two-model waveform-fit comparison plates (matplotlib + obspy only).

Run inside a fwd_* directory after generating waveforms_A/ and waveforms_B/
(each a copy of waveforms/ from an ./ffwd run of one model). Overlays
obs (black) with syn A (crimson) and syn B (steelblue), annotates per-trace
zero-lag correlation for both, and prints a per-band summary table
(mean CC and variance reduction). Writes fits_cmp_{P,Pvel,SH,rayl}.pdf.

usage: plot_fit_comparison_mpl.py [dirA dirB labelA labelB]
       (defaults: waveforms_mcmc waveforms_eks MCMC EKS)
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from obspy import read

PER_COL, NCOL = 12, 3   # 36/page: tighter plates, fewer pages (user 2026-08-13)
LW_OBS, LW_SYN = 2.0, 2.0
COL_A, COL_B = "crimson", "steelblue"

def par(key, path="Par.file"):
    for line in open(path):
        if line.split("=")[0].strip() == key.rstrip("="):
            return float(line.split("=")[1].split("#")[0])
    raise KeyError(key)

def par_or(key, default):
    try:
        return par(key)
    except KeyError:
        return default

def stations(path):
    out = []
    for line in open(path):
        t = line.split()
        if t:
            out.append((".".join(os.path.basename(t[0]).split(".")[:2]),
                        float(t[1]), float(t[2])))
    return out

def cc(o, s):
    n = min(len(o), len(s))
    o, s = o[:n], s[:n]
    d = np.sqrt(np.dot(o, o) * np.dot(s, s))
    return float(np.dot(o, s) / d) if d > 0 else 0.0

class BandStats:
    def __init__(self):
        self.ccs, self.num, self.den = [], 0.0, 0.0
    def add(self, o, s):
        n = min(len(o), len(s))
        o, s = o[:n], s[:n]
        self.ccs.append(cc(o, s))
        self.num += float(np.dot(o - s, o - s))
        self.den += float(np.dot(o, o))
    def row(self):
        return (np.mean(self.ccs) if self.ccs else np.nan,
                1.0 - self.num / self.den if self.den > 0 else np.nan)

STATS = {}  # (band, model) -> BandStats

def plate(pdfname, entries, bg, nd, title, band, wa, wb, la, lb, cut=0.0):
    """entries: list of (label, sublabel, obs, synA, synB) sorted for display"""
    if not entries:
        return
    # BALANCED pagination: never a nearly-empty trailing page (56 traces ->
    # 19/19/18, not 27/27/2); per-column count shrinks accordingly
    per_page = PER_COL * NCOL
    npages = max(1, -(-len(entries) // per_page))
    per_page = -(-len(entries) // npages)
    pages = [entries[i:i + per_page] for i in range(0, len(entries), per_page)]
    with PdfPages(pdfname) as pdf:
        for ip, chunk in enumerate(pages):
            fig, axes = plt.subplots(1, NCOL, figsize=(8.5, 11))
            axes = np.atleast_1d(axes)
            for ic in range(NCOL):
                ax = axes[ic]
                percol = -(-per_page // NCOL)
                sub = chunk[ic * percol:(ic + 1) * percol]
                for k, (lab, sub2, fo, fa, fb) in enumerate(sub):
                    y0 = len(sub) - k
                    try:
                        o = read(fo)[0]; sa = read(fa)[0]; sb = read(fb)[0]
                    except Exception:
                        continue
                    t = bg + np.arange(o.stats.npts) * o.stats.delta
                    vis = (t >= bg + cut) & (t <= nd - cut)
                    ca, cb = cc(o.data[vis], sa.data[:len(t)][vis]), \
                             cc(o.data[vis], sb.data[:len(t)][vis])
                    a = max(np.max(np.abs(o.data[vis]), initial=0),
                            np.max(np.abs(sa.data[:len(t)][vis]), initial=0),
                            np.max(np.abs(sb.data[:len(t)][vis]), initial=0), 1e-30)
                    ax.plot(t, o.data / a * 0.42 + y0, color="0.1", lw=LW_OBS)
                    ax.plot(bg + np.arange(sa.stats.npts) * sa.stats.delta,
                            sa.data / a * 0.42 + y0, color=COL_A, lw=LW_SYN)
                    ax.plot(bg + np.arange(sb.stats.npts) * sb.stats.delta,
                            sb.data / a * 0.42 + y0, color=COL_B, lw=LW_SYN)
                    ax.text(bg + cut, y0 + 0.27, lab, fontsize=9.5,
                            va="bottom", fontweight="bold")
                    ax.text(nd - cut, y0 + 0.27,
                            f"{ca:.2f}/{cb:.2f}", fontsize=8.5, va="bottom",
                            ha="right", color="0.3")
                    ax.text(bg + cut, y0 - 0.42, sub2, fontsize=8.5,
                            va="bottom", color="0.25")
                ax.set_xlim(bg + cut, nd - cut)
                ax.set_ylim(0.2, -(-per_page // NCOL) + 1)
                ax.set_yticks([])
                ax.set_xlabel("time (s)", fontsize=10)
                for sp in ("top", "right", "left"):
                    ax.spines[sp].set_visible(False)
                ax.tick_params(labelsize=9)
            fig.suptitle(f"{title} \N{EM DASH} obs (black), {la} (red), {lb} (blue);"
                         f" CC {la}/{lb} at right"
                         + (f"  ({ip+1}/{len(pages)})" if len(pages) > 1 else ""),
                         fontsize=10.5)
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(fig)
            plt.close(fig)
    # band stats over the full (visible) set
    for lab, sub2, fo, fa, fb in entries:
        try:
            o = read(fo)[0]; sa = read(fa)[0]; sb = read(fb)[0]
        except Exception:
            continue
        t = bg + np.arange(o.stats.npts) * o.stats.delta
        vis = (t >= bg + cut) & (t <= nd - cut)
        STATS.setdefault((band, "A"), BandStats()).add(o.data[vis], sa.data[:len(t)][vis])
        STATS.setdefault((band, "B"), BandStats()).add(o.data[vis], sb.data[:len(t)][vis])
    print(f"wrote {pdfname} ({len(entries)} traces, {len(pages)} page(s))")

def main():
    wa = sys.argv[1] if len(sys.argv) > 1 else "waveforms_mcmc"
    wb = sys.argv[2] if len(sys.argv) > 2 else "waveforms_eks"
    la = sys.argv[3] if len(sys.argv) > 3 else "MCMC"
    lb = sys.argv[4] if len(sys.argv) > 4 else "EKS"
    cutP = par_or("maxshftP", 2.0)
    cutSH = par_or("maxshftSH", 5.0)
    cutR = par_or("maxshftRayl", 3.0)
    staP = stations("stations.info")
    ndisp = sum(1 for n, _, _ in staP if not n.startswith("vel_"))
    bg, nd = par("bg_timeP"), par("nd_timeP")
    def pent(i, n, d, az):
        return (n.replace("vel_", ""), f"{d:.1f}\N{DEGREE SIGN} az {az:.0f}\N{DEGREE SIGN}",
                f"{wa}/P_obs_{i:04d}.sac", f"{wa}/P_syn_{i:04d}.sac",
                f"{wb}/P_syn_{i:04d}.sac")
    disp = sorted((az, pent(i, n, d, az)) for i, (n, d, az) in enumerate(staP) if i < ndisp)
    vel = sorted((az, pent(i, n, d, az)) for i, (n, d, az) in enumerate(staP) if i >= ndisp)
    plate("fits_cmp_P.pdf", [p for _, p in disp], bg, nd,
          "Teleseismic P (displacement)", "P disp", wa, wb, la, lb, cut=cutP)
    plate("fits_cmp_Pvel.pdf", [p for _, p in vel], bg, nd,
          "Teleseismic P (velocity)", "P vel", wa, wb, la, lb, cut=cutP)
    staSH = stations("stationsSH.info")
    bg, nd = par("bg_timeSH"), par("nd_timeSH")
    sh = sorted((az, (n, f"{d:.1f}\N{DEGREE SIGN} az {az:.0f}\N{DEGREE SIGN}",
                      f"{wa}/SH_obs_{i:04d}.sac", f"{wa}/SH_syn_{i:04d}.sac",
                      f"{wb}/SH_syn_{i:04d}.sac"))
                for i, (n, d, az) in enumerate(staSH))
    plate("fits_cmp_SH.pdf", [p for _, p in sh], bg, nd,
          "Teleseismic SH", "SH", wa, wb, la, lb, cut=cutSH)
    staL = stations("stationsloc.info")
    evlo, evla = par("evlo"), par("evla")
    bg, nd = par("bg_timeRayl"), par("nd_timeRayl")
    def dist_km(lo, ll):
        return np.hypot((lo - evlo) * 111.32 * np.cos(np.radians(evla)),
                        (ll - evla) * 110.574)
    rows = sorted((dist_km(lo, ll), i, n) for i, (n, lo, ll) in enumerate(staL))
    ents = []
    for dk, i, n in rows:
        for comp in "enz":
            ents.append((f"{n}.{comp.upper()}", f"{dk:.0f} km",
                         f"{wa}/rayl{comp}_obs_{i:04d}.sac",
                         f"{wa}/rayl{comp}_syn_{i:04d}.sac",
                         f"{wb}/rayl{comp}_syn_{i:04d}.sac"))
    plate("fits_cmp_rayl.pdf", ents, bg, nd,
          "Regional full waveforms", "Regional", wa, wb, la, lb, cut=cutR)
    print(f"\n{'band':>10} {'n':>4}   mean CC {la}/{lb}    VR {la}/{lb}")
    for band in ("P disp", "P vel", "SH", "Regional"):
        if (band, "A") not in STATS:
            continue
        ca, va = STATS[(band, "A")].row()
        cbb, vb = STATS[(band, "B")].row()
        n = len(STATS[(band, "A")].ccs)
        print(f"{band:>10} {n:>4}   {ca:.3f} / {cbb:.3f}      {va:.3f} / {vb:.3f}")

if __name__ == "__main__":
    main()
