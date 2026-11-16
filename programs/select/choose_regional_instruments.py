#!/usr/bin/env python3
"""Choose one instrument per regional station, and flag clipped broadband records (step 1).

Called by MyExtractSeed_loc.sh inside IRISloc/<event>_loc/:
    python choose_regional_instruments.py clip-check   on the raw records (counts), before processing
    python choose_regional_instruments.py choose       on Vel/enz_<family>/ after every family was processed

Why. Until v3 the regional step used strong-motion HN? channels only, and broadband channels only
when no HN? channel existed at all; HH? was never read. A download with any strong-motion record
lost every broadband-only station: the 2026 Colombia M7.4 kept no station closer than 496 km while
HEL (169 km) and RUS (369 km) had HH? data (cohort report 2026-10-23, item 3); the 2026 Panama M7.7
clean-environment test kept 4 of ~20 stations. Several sensors at one site (IU.MAJO.00, .10, .60
for the 2025 Aomori M7.6) were all kept and counted the site two or three times.

Rule, one record per NET.STA (thresholds measured by the cohort report on the Colombian network):
  * clipped (raw counts, any component): a plateau of >= 3 consecutive samples AT the peak value
    with the peak above 10% of 2^23 counts (a digitizer clip is a flat top at its largest code).
    Two rules were tried and dropped: "a plateau above half the peak" flagged 10 of 22 unclipped
    strong-motion stations of the 2026 Venezuela example (smooth, oversampled 200-sps records repeat
    values near a peak: CM.ARGC HNE 6 equal samples at 0.6% of 2^23), and "peak >= 95% of the 24-bit
    full scale" assumes 24-bit digitizers: IU.PET in the Kamchatka example peaks at 4.8-6.6 x 2^23
    counts without any flat top. A saturated SENSOR (not digitizer) is caught by the agreement test
    below: CM.OCA.00.HH disagreed with its strong-motion twin (cc 0-0.5, amplitude x2.5 at
    0.02-0.05 Hz) while every other broadband record matched its twin (cc >= 0.92, ratio 0.9-1.1).
  * broadband (HH before BH, lowest location code first) and strong motion (HN) both present: keep
    the unclipped broadband record if it agrees with the strong-motion one on every component
    (cc >= 0.9 and peak ratio 0.8-1.25 at 0.02-0.1 Hz, the regional inversion band); otherwise keep
    the strong-motion record. Broadband is preferred because accelerometers carry tilt and baseline
    drift at the longest periods of the band.
  * one instrument type only: keep it unless it is clipped (a clipped record is dropped: a zero
    component weight would leave a station without one component, which the weighting does not
    handle).
The other records move to Vel/not_selected/; decisions are written to regional_instrument_choice.log.
"""
import csv
import glob
import os
import shutil
import sys

import numpy as np
from obspy import read

flat_run_minimum = 3
plateau_minimum_fraction = 0.1
full_scale_counts = 2 ** 23
correlation_minimum = 0.9
amplitude_ratio_range = (0.8, 1.25)
agreement_band_hz = (0.02, 0.1)
clip_table = "regional_clip_check.csv"
component_of_channel = {"Z": "z", "E": "e", "N": "n", "1": "e", "2": "n"}   # 1/2 are rotated to e/n


def longest_flat_run_at_peak(samples):
    """length of the longest run of consecutive samples equal to +peak or -peak (1 for an ordinary peak)"""
    samples = np.asarray(samples, dtype=np.float64)
    if samples.size < 2:
        return 0
    peak = np.max(np.abs(samples))
    same_as_previous = np.concatenate([[False], samples[1:] == samples[:-1]]) & (np.abs(samples) == peak)
    longest = current = 0
    for is_same in same_as_previous:
        current = current + 1 if is_same else 0
        longest = max(longest, current)
    return longest + 1 if longest else 0


def clip_check():
    """raw records (counts) in the current folder, named NET.STA.LOC.CHA by MyExtractSeed_loc.sh"""
    rows = []
    for file_name in sorted(glob.glob("*.*.*.[HB][HN]?")):
        network, station, location, channel = file_name.split(".")[:4]
        try:
            samples = read(file_name)[0].data
        except Exception as error:
            print(f"choose_regional_instruments: cannot read {file_name}: {error}")
            continue
        peak = float(np.max(np.abs(samples))) if samples.size else 0.0
        flat_run = longest_flat_run_at_peak(samples)
        clipped = flat_run >= flat_run_minimum and peak >= plateau_minimum_fraction * full_scale_counts
        rows.append(dict(network=network, station=station, location=location, family=channel[:2],
                         component=component_of_channel.get(channel[2], channel[2].lower()), peak_counts=int(peak),
                         fraction_of_full_scale=round(peak / full_scale_counts, 3), flat_run=flat_run, clipped=int(clipped)))
    with open(clip_table, "w", newline="") as table:
        writer = csv.DictWriter(table, fieldnames=list(rows[0]) if rows else ["network"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"choose_regional_instruments: clip check of {len(rows)} raw components, "
          f"{sum(row['clipped'] for row in rows)} clipped -> {clip_table}")


def band_limited(file_name):
    trace = read(file_name)[0]
    trace.detrend("demean")
    trace.detrend("linear")
    trace.taper(0.05)
    trace.filter("bandpass", freqmin=agreement_band_hz[0], freqmax=agreement_band_hz[1], corners=4, zerophase=True)
    return trace


def agreement(broadband_stem, strong_motion_stem):
    """per component: (correlation, peak ratio broadband/strong motion) on the common time window"""
    out = {}
    for component in "zne":
        broadband = band_limited(f"{broadband_stem}.{component}")
        strong_motion = band_limited(f"{strong_motion_stem}.{component}")
        start = max(broadband.stats.starttime, strong_motion.stats.starttime) + 20
        end = min(broadband.stats.endtime, strong_motion.stats.endtime) - 20
        broadband.trim(start, end)
        strong_motion.trim(start, end)
        npts = min(broadband.stats.npts, strong_motion.stats.npts)
        x, y = broadband.data[:npts], strong_motion.data[:npts]
        energy = np.sqrt(np.dot(x, x) * np.dot(y, y))
        correlation = float(np.dot(x, y) / energy) if energy > 0 else 0.0
        ratio = float(np.max(np.abs(x)) / np.max(np.abs(y))) if np.max(np.abs(y)) > 0 else 0.0
        out[component] = (correlation, ratio)
    return out


def choose():
    clipped = {}
    if os.path.exists(clip_table):
        for row in csv.DictReader(open(clip_table)):
            key = (row["network"], row["station"], row["location"], row["family"])
            clipped.setdefault(key, set())
            if int(row["clipped"]):
                clipped[key].add(row["component"])
    candidates = {}
    for z_file in glob.glob("Vel/enz_*/*.*.*.z"):
        family = os.path.basename(os.path.dirname(z_file))[len("enz_"):]
        stem = z_file[:-2]
        if not all(os.path.exists(f"{stem}.{component}") for component in "en"):
            continue
        network, station, location = os.path.basename(stem).split(".")[:3]
        candidates.setdefault((network, station), []).append(
            dict(stem=stem, family=family, location=location,
                 clipped="".join(sorted(clipped.get((network, station, location, family), set())))))
    os.makedirs("Vel/enz_selected", exist_ok=True)
    os.makedirs("Vel/not_selected", exist_ok=True)
    decisions = []
    for (network, station), records in sorted(candidates.items()):
        broadband = sorted((r for r in records if r["family"] in ("HH", "BH")), key=lambda r: (r["family"] != "HH", r["location"]))
        strong_motion = sorted((r for r in records if r["family"] == "HN"), key=lambda r: r["location"])
        unclipped_broadband = [r for r in broadband if not r["clipped"]]
        unclipped_strong_motion = [r for r in strong_motion if not r["clipped"]]
        chosen, reason = None, ""
        if unclipped_broadband and strong_motion:
            measured = agreement(unclipped_broadband[0]["stem"], strong_motion[0]["stem"])
            summary = " ".join(f"{c}:cc{cc:.2f}/x{ratio:.2f}" for c, (cc, ratio) in measured.items())
            agrees = all(cc >= correlation_minimum and amplitude_ratio_range[0] <= ratio <= amplitude_ratio_range[1]
                         for cc, ratio in measured.values())
            if agrees or not unclipped_strong_motion:
                chosen, reason = unclipped_broadband[0], f"broadband agrees with strong motion ({summary})" if agrees else \
                    f"strong motion clipped; broadband kept although they disagree ({summary})"
            else:
                chosen, reason = unclipped_strong_motion[0], f"broadband disagrees with strong motion ({summary})"
        elif unclipped_broadband:
            chosen, reason = unclipped_broadband[0], "broadband only"
        elif unclipped_strong_motion:
            chosen, reason = unclipped_strong_motion[0], ("broadband clipped" if broadband else "strong motion only")
        for record in records:
            target = "Vel/enz_selected" if record is chosen else "Vel/not_selected"
            for component in "enz":
                name = os.path.basename(f"{record['stem']}.{component}")
                if record is not chosen:
                    name = f"{record['family']}_{name}"
                shutil.move(f"{record['stem']}.{component}", os.path.join(target, name))
        if chosen is None:
            decisions.append(f"{network}.{station}: DROPPED - every record clipped "
                             + ", ".join(f"{r['family']}.{r['location']} ({r['clipped']})" for r in records))
        else:
            decisions.append(f"{network}.{station}: {chosen['family']} location '{chosen['location']}' - {reason}"
                             + (f" (of {len(records)} records)" if len(records) > 1 else ""))
    open("regional_instrument_choice.log", "w").write("\n".join(decisions) + "\n")
    print("\n".join(decisions))


if __name__ == "__main__":
    {"clip-check": clip_check, "choose": choose}[sys.argv[1]]()
