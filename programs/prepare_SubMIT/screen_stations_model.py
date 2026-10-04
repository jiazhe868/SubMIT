#!/usr/bin/env python3
"""Stage-2 MODEL-BASED station screen (user failure-mode taxonomy 2026-08-11).

Run from an IRIS dir after the 1-sub scan produced a full (non-compact) ffwd
run in fwd_<event>/ (sta_residual.dat + waveforms/). Two screens that no
data-only rule can implement:

- POLARITY FLIP (2024 Chile: GT.DBIC, CN.DRLN, GE.SNAA): the per-trace
  best-alignment CC from sta_residual.dat is strongly NEGATIVE - the trace
  anticorrelates with the prediction no matter the shift.
- GAIN/RESPONSE ERROR (2024 Chile: CN.GAC): peak |obs|/|syn| deviates from
  the band's median ratio by more than GAIN_FACTOR - shape fine, scale wrong.

The 1-sub model is crude but these diagnostics are robust to it: polarity and
order-of-magnitude gain do not depend on source complexity. Writes
excluded_stations_model.txt and moves flagged traces to data/excluded (and
dataloc/excluded) in EVERY inv_<event>_*sub dir; rerun prepare_stationinfo
and the weight balancing afterwards (step3 hook does both).

usage: screen_stations_model.py <event_name>
"""
import glob
import os
import shutil
import sys

import numpy as np
from obspy import read

CC_FLIP = 0.0        # best-alignment CC that is NEGATIVE at every shift =
                     # polarity/broken (2024 Chile: CN.DRLN -0.20, GE.SNAA
                     # -0.05; worst legitimate station +0.16)
GAIN_FACTOR = 4.0    # |log10 obs/syn| beyond log10(4) from band median

ev = sys.argv[1]
fwd = f"fwd_{ev}"

def parse_sta_residual(path):
    """returns {band: [(idx, resid, cc, tshift), ...]}"""
    bands, cur = {}, None
    for line in open(path):
        t = line.split()
        if line.startswith("#"):
            cur = line.strip("# \n")
            bands[cur] = []
        elif len(t) >= 4 and cur:
            bands[cur].append((int(t[0]), float(t[1]), float(t[2]), float(t[3])))
    return bands

def station_files(info):
    out = []
    for line in open(os.path.join(fwd, info)):
        t = line.split()
        if t:
            out.append(t[0])          # e.g. data/IU.ANMO.00.z
    return out

def peak(f):
    try:
        return float(np.max(np.abs(read(f)[0].data)))
    except Exception:
        return 0.0

resid = parse_sta_residual(os.path.join(fwd, "sta_residual.dat"))
report, flagged = [], set()

for band, info, obspat, synpat, comp in (
        ("P waves", "stations.info", "P_obs_%04d.sac", "P_syn_%04d.sac", "z"),
        ("SH waves", "stationsSH.info", "SH_obs_%04d.sac", "SH_syn_%04d.sac", "t")):
    files = station_files(info)
    rows = resid.get(band, [])
    ratios = {}
    for idx, _, cc, _ in rows:
        j = idx - 1
        if j >= len(files):
            continue
        f = files[j]
        if cc < CC_FLIP:
            report.append(f"MODEL {band} {f} best CC {cc:.2f} < {CC_FLIP} (polarity?)")
            flagged.add(f)
        o = peak(os.path.join(fwd, "waveforms", obspat % j))
        sy = peak(os.path.join(fwd, "waveforms", synpat % j))
        if o > 0 and sy > 0:
            ratios[f] = np.log10(o / sy)
    if ratios and not os.environ.get('STAGE2_POLARITY_ONLY'):
        med = float(np.median(list(ratios.values())))
        for f, r in ratios.items():
            if abs(r - med) > np.log10(GAIN_FACTOR):
                report.append(f"MODEL {band} {f} obs/syn x{10**(r-med):.1f} of band "
                              f"median (gain/response?)")
                flagged.add(f)

if not flagged:
    print("screen_stations_model: nothing flagged")
    sys.exit(0)

# expand: P disp <-> vel twins share fate; write list; apply to every inv dir
expanded = set(flagged)
for f in list(flagged):
    base = os.path.basename(f)
    twin = "data/" + (base[4:] if base.startswith("vel_") else "vel_" + base)
    if os.path.exists(os.path.join(fwd, twin)):
        expanded.add(twin)

with open("excluded_stations_model.txt", "w") as fh:
    fh.write("\n".join(sorted(report)) + "\n")
for line in report:
    print(line)

napplied = 0
for d in sorted(glob.glob(f"inv_{ev}_*sub")) + [fwd]:
    for f in sorted(expanded):
        src = os.path.join(d, f)
        if os.path.exists(src):
            dst = os.path.join(d, os.path.dirname(f), "excluded")
            os.makedirs(dst, exist_ok=True)
            shutil.move(src, os.path.join(dst, os.path.basename(f)))
            napplied += 1
print(f"screen_stations_model: {len(expanded)} trace(s) excluded across dirs "
      f"({napplied} moves)")
