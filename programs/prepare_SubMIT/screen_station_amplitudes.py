"""Per-station amplitude screening for SubMIT (run inside an inv_* folder, after
process_tel_comp.sh and BEFORE prepare_stationinfo.sh / gen_weight.py).

SNR-based selection upstream cannot catch stations whose waveforms are clean but
whose amplitude is wrong (bad instrument response, array elements) or near-nodal:
e.g. 2026 Calama M6.9, IM.PD31 3.4x larger and AU.MAW 3x smaller than any other
station, both dominating the L2 misfit. Here each wave-type group (P displacement
*.z, P velocity vel_*.z, SH *.t) is screened on log10 peak amplitude with a
robust median/MAD rule; outliers move to data/excluded/ (never deleted) and are
logged in excluded_stations.txt.
"""
import glob
import os
import shutil
import numpy as np
from obspy import read

MAD_K = 3.0          # flag beyond median +/- MAD_K * 1.4826 * MAD ...
MIN_FACTOR = 8.0     # ... but never flag within a factor MIN_FACTOR of the median
MIN_GROUP = 8        # do not screen groups with fewer stations than this
MIN_GROUP_LOC = 4    # regional sets are small; still screen (factor rule catches
                     # gross unit errors like raw counts, e.g. C1.TA01 ~1e5 too big)

def peak_log_amp(f):
    tr = read(f)[0]
    a = np.max(np.abs(tr.data))
    return np.log10(a) if a > 0 else -np.inf

def screen(files, label, report, min_group=MIN_GROUP):
    files = sorted(files)
    if len(files) < min_group:
        return []
    amps = np.array([peak_log_amp(f) for f in files])
    ok = np.isfinite(amps)
    med = np.median(amps[ok])
    mad = np.median(np.abs(amps[ok] - med))
    thr = max(MAD_K * 1.4826 * mad, np.log10(MIN_FACTOR))
    out = []
    for f, a in zip(files, amps):
        if not np.isfinite(a) or abs(a - med) > thr:
            ratio = 10 ** (a - med) if np.isfinite(a) else 0.0
            report.append(f"{label} {f} log10amp={a:.2f} median={med:.2f} "
                          f"(x{ratio:.2f} of median)")
            out.append(f)
    return out

SNR_MIN_TEL = 5.5    # SNR floor for genuinely NOISE-DOMINATED tel traces
                     # (2024 California II.KWJN P at 5.1 also junk)
                     # only (2024 Chile: G.TAOE 1.1/1.3, II.ASCN 3.8). Gain,
                     # polarity and drift failures have NORMAL SNR and are
                     # handled by dedicated screens (drift below; polarity &
                     # gain by the model-based stage-2 screen after the 1-sub
                     # scan) - a high SNR floor mislabels those modes.
SNR_MIN_LOC = 3.0    # regional floor
DRIFT_LF_FRAC = 0.4  # REGIONAL drift screen: fraction of spectral energy
                     # below 0.04 Hz (2024 Chile: C1.A03C.n 0.96, C1.A01C.z
                     # 0.74; next healthy station 0.21). Regional only -
                     # teleseismic displacement legitimately has long-period
                     # energy.

def drift_screen(files, label, report):
    """ultra-long-period drift: one giant >25-s swing dominating a regional
    velocity record (instrument tilt/response drift)"""
    out = []
    try:
        from obspy import read
    except Exception:
        return out
    for f in files:
        try:
            tr = read(f)[0]
            d = tr.data.astype(float)
            d -= d.mean()
            if np.max(np.abs(d)) <= 0:
                continue
            F = np.abs(np.fft.rfft(d * np.hanning(len(d)))) ** 2
            fr = np.fft.rfftfreq(len(d), tr.stats.delta)
            lf = float(F[fr < 0.04].sum() / max(F.sum(), 1e-30))
            if lf > DRIFT_LF_FRAC:
                report.append(f"DRIFT {label} {f} lowfreq_frac={lf:.2f} > {DRIFT_LF_FRAC}")
                out.append(f)
        except Exception:
            continue
    return out

def snr_screen(files, label, report, snr_min, sig_end=90):
    """signal-to-noise + dead-channel screen: noise = RMS of the pre-arrival
    window (t1-120..t1-10 s, clipped to trace start), signal = peak |.| in
    t1-5..t1+90 s. A channel with EXACTLY zero pre-arrival noise is a broken
    record (2024 Chile C1.A03C: zero-padded start + 18x inter-component
    amplitude inconsistency) and is rejected outright."""
    out = []
    try:
        from obspy import read
    except Exception:
        return out
    for f in files:
        try:
            tr = read(f)[0]
            b = tr.stats.sac.b
            t1 = tr.stats.sac.get("t1", None)
            dt = tr.stats.delta
            if t1 is None:
                continue
            i1 = int(round((t1 - b) / dt))
            noise = tr.data[max(0, i1 - int(120/dt)):max(1, i1 - int(10/dt))].astype(float)
            # sig_end=None: signal window runs to END of trace - REGIONAL
            # surface waves arrive at ~dist/3 s (2026 Venezuela: 100-230 s for
            # the 300-700 km Colombian stations; a fixed +90 s window measured
            # pre-surface-wave quiet as "signal" and mass-excluded the band)
            iend = tr.stats.npts if sig_end is None else i1 + int(sig_end/dt)
            sig = tr.data[max(0, i1 - int(5/dt)):iend].astype(float)
            if len(noise) < 10 or len(sig) < 10:
                continue
            nr = float(np.sqrt(np.mean(noise ** 2)))
            if nr == 0.0:
                report.append(f"SNR {label} {f} DEAD pre-arrival window (broken record)")
                out.append(f)
                continue
            snr = float(np.max(np.abs(sig))) / nr
            if snr < snr_min:
                report.append(f"SNR {label} {f} snr={snr:.1f} < {snr_min}")
                out.append(f)
        except Exception:
            continue
    return out

def main():
    os.makedirs("data/excluded", exist_ok=True)
    report = []
    zdisp = [f for f in glob.glob("data/*.z")
             if not os.path.basename(f).startswith("vel_")]
    zvel = glob.glob("data/vel_*.z")
    tsh = glob.glob("data/*.t")

    bad = set(screen(zdisp, "P-disp", report))
    bad |= set(screen(zvel, "P-vel", report))
    bad |= set(screen(tsh, "SH", report))
    bad |= set(snr_screen(zdisp, "P-disp", report, SNR_MIN_TEL))
    bad |= set(snr_screen(tsh, "SH", report, SNR_MIN_TEL))

    # regional 3-component sets: screen all components together; exclude whole station
    loc = sorted(glob.glob("dataloc/*.e") + glob.glob("dataloc/*.n") + glob.glob("dataloc/*.z"))
    badloc = set(screen(loc, "LOC", report, min_group=MIN_GROUP_LOC))
    badloc |= set(snr_screen(loc, "LOC", report, SNR_MIN_LOC, sig_end=None))
    badloc |= set(drift_screen(loc, "LOC", report))
    for f in list(badloc):
        stem = f[:-1]              # strip component letter
        for c in "enz":
            if os.path.exists(stem + c) and stem + c not in badloc:
                report.append(f"pair    {stem + c} (component of excluded station)")
                badloc.add(stem + c)
    if badloc:
        os.makedirs("dataloc/excluded", exist_ok=True)
        for f in sorted(badloc):
            shutil.move(f, os.path.join("dataloc/excluded", os.path.basename(f)))

    # keep displacement/velocity P pairs consistent: exclude both if either fails
    for f in list(bad):
        base = os.path.basename(f)
        twin = ("data/" + base[4:]) if base.startswith("vel_") else ("data/vel_" + base)
        if os.path.exists(twin) and twin not in bad:
            report.append(f"pair    {twin} (companion of excluded {f})")
            bad.add(twin)

    for f in sorted(bad):
        shutil.move(f, os.path.join("data/excluded", os.path.basename(f)))

    # expert QC that no automatic rule captures (near-nodal geometry, path
    # complexity, visually bad fits): list SAC basenames - one per line, "#"
    # comments allowed - in excluded_stations_manual.txt (tel names move from
    # data/, names containing "loc:" prefix move from dataloc/)
    if os.path.exists("excluded_stations_manual.txt"):
        nman = 0
        for line in open("excluded_stations_manual.txt"):
            name = line.split("#")[0].strip()
            if not name:
                continue
            if name.startswith("loc:"):
                src, dstd = os.path.join("dataloc", name[4:]), "dataloc/excluded"
            else:
                src, dstd = os.path.join("data", name), "data/excluded"
            hits = glob.glob(src) if any(c in name for c in "*?") else \
                ([src] if os.path.exists(src) else [])
            for h in hits:
                os.makedirs(dstd, exist_ok=True)
                shutil.move(h, os.path.join(dstd, os.path.basename(h)))
                report.append(f"manual  {h}")
                nman += 1
        if nman:
            print(f"screen_station_amplitudes: {nman} file(s) excluded by "
                  f"excluded_stations_manual.txt")

    # clipped regional recordings (giant events saturate nearby broadbands):
    # a clipped trace has a plateau at its extremes - flag when >1% of samples
    # sit within 1% of the absolute peak
    for f in sorted(glob.glob("dataloc/*.[enz]")):
        try:
            tr = read(f)[0]
            a = np.abs(tr.data.astype(float))
            peak = a.max()
            if peak <= 0:
                continue
            frac = float(np.mean(a > 0.99 * peak))
            if frac > 0.01:
                stem = f[:-1]
                os.makedirs("dataloc/excluded", exist_ok=True)
                for c in "enz":
                    if os.path.exists(stem + c):
                        shutil.move(stem + c, os.path.join("dataloc/excluded",
                                    os.path.basename(stem + c)))
                        report.append(f"clipped {stem + c} (plateau fraction {frac:.3f})")
        except Exception:
            pass

    with open("excluded_stations.txt", "w") as fp:
        fp.write("\n".join(report) + ("\n" if report else ""))
    print(f"screen_station_amplitudes: excluded {len(bad)} tel + {len(badloc)} loc "
          f"file(s) of {len(zdisp) + len(zvel) + len(tsh) + len(loc)}; "
          f"see excluded_stations.txt")

if __name__ == "__main__":
    main()
