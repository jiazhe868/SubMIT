#!/usr/bin/env python3
"""Write the station lists that forward_model.c reads (run inside an inversion folder).

  stations.info     data/<record>.z then data/vel_<record>.z: distance (deg), azimuth (deg),
                    P and S apparent velocities 111.195/ray parameter (km/s), and the GF distance key
                    (integer km; vel_<km> for the velocity records)
  stationsSH.info   the same for data/<record>.t
  stationsloc.info  dataloc/<record>. (components e, n, z appended by the reader), longitude,
                    latitude, NET.STA, and the weights 1 1 1 of the e, n, z components; only
                    stations within nd_localdist - 200 km (the regional GF grid keeps a 200-km margin
                    for subevents away from the epicentre), sorted by distance

Why Python: v2 launched the Java TauP twice per station (prepare_stationinfo.sh +
doslownesstaup.sh + slownesstaup.sh + slownesstaupS.sh), about 900 launches for a 220-station
event, and the step-3 calibration reruns this in every inversion folder. ObsPy's TauP port reads
the same iasp91 model; the first P (or p) and S (or s) arrival give the ray parameters as before.

usage: python prepare_stationinfo.py
"""
import glob
import os

from obspy import read
from obspy.taup import TauPyModel

kilometres_per_degree = 111.195
model = TauPyModel("iasp91")
ray_parameter_cache = {}


def apparent_velocity(depth_km, distance_degrees, phase):
    """111.195 / ray parameter (s/deg) of the first arrival of phase (P/p or S/s), km/s"""
    key = (round(depth_km, 3), round(distance_degrees, 4), phase)
    if key not in ray_parameter_cache:
        names = [phase, phase.lower()]
        arrivals = model.get_travel_times(source_depth_in_km=max(depth_km, 0.0), distance_in_degree=distance_degrees,
                                          phase_list=names)
        arrivals = [arrival for arrival in arrivals if arrival.name in names]
        if not arrivals:
            raise SystemExit(f"prepare_stationinfo: no {phase} arrival at {distance_degrees:.2f} deg")
        ray_parameter_cache[key] = kilometres_per_degree / arrivals[0].ray_param_sec_degree
    return ray_parameter_cache[key]


def green_function_distance_key(distance_km):
    """the distance part of the Green's-function file names, computed exactly as gf_SubMIT/dotel.sh
    computes it: saclst prints dist with 7 significant digits and gawk's %d truncates that print.
    int(distance_km) alone differs whenever dist lies just below an integer (6202.9995 -> 6202, while
    dotel.sh wrote 6203.grn.*: 3 of 309 test values), and step 3 then stops on a missing file."""
    return int(float(f"{distance_km:.7g}"))


def teleseismic_rows(file_names):
    """one row per record: name, distance, azimuth, c_p, c_s, GF distance key"""
    rows = []
    for file_name in file_names:
        header = read(file_name, headonly=True)[0].stats.sac
        rows.append((file_name, header.gcarc, header.az, apparent_velocity(header.evdp, header.gcarc, "P"),
                     apparent_velocity(header.evdp, header.gcarc, "S"), green_function_distance_key(header.dist)))
    return rows


def write_rows(path, rows):
    with open(path, "w") as out:
        for name, distance, azimuth, p_velocity, s_velocity, key in rows:
            out.write(f"{name} {distance:.5f} {azimuth:.5f} {p_velocity:.4f} {s_velocity:.4f} {key}\n")


displacement_records = sorted(f for f in glob.glob("data/*.z") if not os.path.basename(f).startswith("vel_"))
p_rows = teleseismic_rows(displacement_records)
velocity_rows = [(name.replace("data/", "data/vel_", 1), distance, azimuth, p_velocity, s_velocity, f"vel_{key}")
                 for name, distance, azimuth, p_velocity, s_velocity, key in p_rows]
write_rows("stations.info", p_rows + velocity_rows)
write_rows("vel_stations.info", velocity_rows)
write_rows("stationsSH.info", teleseismic_rows(sorted(glob.glob("data/*.t"))))

# regional: the GF grid ends at nd_localdist (distdep.dat column 2), 200 km of which are margin
regional_cutoff_km = 650 - 200
if os.path.exists("distdep.dat"):
    regional_cutoff_km = float(open("distdep.dat").read().split()[1]) - 200
regional = []
for file_name in sorted(glob.glob("dataloc/*.z")):
    header = read(file_name, headonly=True)[0].stats.sac
    if header.dist < regional_cutoff_km:
        network_station = ".".join(os.path.basename(file_name).split(".")[:2])
        regional.append((header.dist, f"{file_name[:-1]} {header.stlo:g} {header.stla:g} {network_station} 1 1 1"))
with open("stationsloc.info", "w") as out:
    for _, line in sorted(regional, key=lambda item: item[0]):
        out.write(line + "\n")
print(f"prepare_stationinfo: {len(p_rows)} P, {len(open('stationsSH.info').readlines())} SH, "
      f"{len(regional)} regional stations (within {regional_cutoff_km:.0f} km)")
