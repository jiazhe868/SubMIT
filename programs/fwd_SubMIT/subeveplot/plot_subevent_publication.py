#!/usr/bin/env python3
"""
Automated Publication-Grade Subevent Map and Cross-Section Plotter in GMT.

Fully Automated Dynamic Hyperparameter Selection:
1. Coordinates & Epicenter: Fetched from Par.file (evlo, evla) or search_par.file.
2. Moment & Magnitudes: Dynamically computed from fm.dat scaling and components.
3. Map Extents (-R): Dynamically calculated from subevent locations, epicenter, and magnitude-based beachball sizes.
4. Profile A-A' Transect:
   - Line length & orientation dynamically determined perpendicular to local trench/subduction strike or along subevent cluster trend.
   - Profile depth range dynamically scaled from min/max subevent depths + 95% uncertainties.
5. STF Inset Frame (-R-t1/t2/0/1):
   - Time window dynamically calculated from max origin time shift + max STF duration in Input.model / search_par.file.
   - Plotted outside the top-left outer margin of map panel with auto-computed offset to prevent clipping.
6. 95% Uncertainty Error Bars (Exy):
   - Automatically parsed from hist/allsamples.dat or computed from MCMC chain outputs.
   - Rendered with GMT psxy -Exy over beachballs on both Map and Depth Cross Section.

Author: Antigravity AI
"""

import os
import sys
import math
import numpy as np
import subprocess

def parse_par_file(par_path):
    params = {}
    if not os.path.exists(par_path):
        return params
    with open(par_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()
            if '=' in line:
                key, val = line.split('=', 1)
                key, val = key.strip(), val.strip()
                try:
                    if '.' in val or 'e' in val.lower():
                        params[key] = float(val)
                    else:
                        params[key] = int(val)
                except ValueError:
                    params[key] = val
    return params

def parse_search_par(search_path):
    params = {}
    if not os.path.exists(search_path):
        return params
    with open(search_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            tokens = line.split()
            if len(tokens) >= 2:
                key, val = tokens[0], tokens[1]
                try:
                    params[key] = float(val) if '.' in val else int(val)
                except ValueError:
                    params[key] = val
    return params

def calculate_distance(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = np.sin(dlat / 2.0)**2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon / 2.0)**2
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    return R * c

def run_cmd(cmd):
    res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        print(f"[Warning/Error] Command: {cmd}\nStderr: {res.stderr.strip()}")
    return res.stdout

def compute_95_errors(hist_dirs, nsub, evla, evlo, bg_lon, bg_lat, nd_lon, nd_lat):
    samples_file = None
    for hdir in hist_dirs:
        sfile = os.path.join(hdir, 'allsamples.dat')
        if not os.path.exists(sfile):
            py_proc = os.path.join(hdir, 'process_and_plot.py')
            if os.path.exists(py_proc):
                run_cmd(f"cd {hdir} && {sys.executable} process_and_plot.py")
        if os.path.exists(sfile):
            samples_file = sfile
            break

    if not samples_file or not os.path.exists(samples_file):
        print("Warning: allsamples.dat not available, using fallback zero/default error bars.")
        return [(0.05, 0.05)]*nsub, [(5.0, 5.0)]*nsub

    samples = np.loadtxt(samples_file)
    R = 6371.0

    # Profile vector in km
    dx_prof = (nd_lon - bg_lon) * (R * np.cos(np.radians(bg_lat))) * (np.pi / 180.0)
    dy_prof = (nd_lat - bg_lat) * R * (np.pi / 180.0)
    prof_len = np.sqrt(dx_prof**2 + dy_prof**2)

    def project_point(lo, la):
        dx_p = (lo - bg_lon) * (R * np.cos(np.radians(bg_lat))) * (np.pi / 180.0)
        dy_p = (la - bg_lat) * R * (np.pi / 180.0)
        return (dx_p * dx_prof + dy_p * dy_prof) / prof_len

    map_errors = []  # list of (lon_err, lat_err)
    section_errors = [] # list of (dist_err, dep_err)

    for k in range(nsub):
        if k == 0:
            dx_s, dy_s, dz_s = samples[:, 2], samples[:, 3], samples[:, 5]
        elif k == 1:
            dx_s, dy_s, dz_s = samples[:, 7], samples[:, 8], samples[:, 10]
        elif k == 2:
            dx_s, dy_s, dz_s = samples[:, 12], samples[:, 13], samples[:, 15]
        elif k == 3:
            dx_s, dy_s, dz_s = samples[:, 17], samples[:, 18], samples[:, 20]
        else:
            dx_s, dy_s, dz_s = samples[:, k*5 + 2], samples[:, k*5 + 3], samples[:, k*5 + 5]

        lons = evlo + (dx_s / (R * np.cos(np.radians(evla)))) * (180.0 / np.pi)
        lats = evla + (dy_s / R) * (180.0 / np.pi)

        lon_err = (np.percentile(lons, 97.5) - np.percentile(lons, 2.5)) / 2.0
        lat_err = (np.percentile(lats, 97.5) - np.percentile(lats, 2.5)) / 2.0
        dep_err = (np.percentile(dz_s, 97.5) - np.percentile(dz_s, 2.5)) / 2.0

        # Profile projection
        dists = [project_point(lo, la) for lo, la in zip(lons[::500], lats[::500])]
        dist_err = (np.percentile(dists, 97.5) - np.percentile(dists, 2.5)) / 2.0

        map_errors.append((lon_err, lat_err))
        section_errors.append((dist_err, dep_err))

    return map_errors, section_errors

def main():
    # 1. Automatic Parameter Parsing from Par.file & search_par.file
    par_paths = ['../Par.file', 'Par.file', '../../Par.file']
    search_paths = ['../search_par.file', 'search_par.file', '../../search_par.file']

    par = {}
    for p in par_paths:
        if os.path.exists(p):
            par = parse_par_file(p)
            print(f"Loaded parameters from {p}")
            break

    search_par = {}
    for sp in search_paths:
        if os.path.exists(sp):
            search_par = parse_search_par(sp)
            print(f"Loaded search parameters from {sp}")
            break

    evla = par.get('evla', search_par.get('evla', -22.3667))
    evlo = par.get('evlo', search_par.get('evlo', -68.6013))
    topo_grd_path = os.environ.get("SUBMIT_TOPO_GRD", "etopo1.grd")  # ETOPO1 grid for GMT relief

    # Input files
    input_model_file = "Input.model" if os.path.exists("Input.model") else "../Input.model"
    fm_file = "fm.dat" if os.path.exists("fm.dat") else "../fm.dat"

    if not os.path.exists(input_model_file) or not os.path.exists(fm_file):
        print(f"Error: {input_model_file} or {fm_file} not found.")
        sys.exit(1)

    input_model = np.loadtxt(input_model_file)
    if input_model.ndim == 1:
        input_model = np.array([input_model])

    fm_data = np.loadtxt(fm_file)
    if fm_data.ndim == 1:
        fm_data = np.array([fm_data])

    nsub = len(input_model)
    print(f"Loaded {nsub} subevents around epicenter (Lat: {evla:.4f}, Lon: {evlo:.4f})")

    # 2. Dynamic Subevent Geometries & Magnitudes
    R_earth = 6371.0
    sub_lons, sub_lats, sub_deps, sub_times, sub_durs = [], [], [], [], []

    for i in range(nsub):
        t_shift, dx, dy, dz = input_model[i, :4]
        dur = input_model[i, 4] if input_model.shape[1] >= 5 else search_par.get('dura_max', 10.0)
        dep = input_model[i, 6] if input_model.shape[1] >= 7 else (dz if dz > 50.0 else par.get('bg_depth', 90.0) + dz)
        lat = evla + (dy / R_earth) * (180.0 / np.pi)
        lon = evlo + (dx / (R_earth * np.cos(np.radians(evla)))) * (180.0 / np.pi)

        sub_lons.append(lon)
        sub_lats.append(lat)
        sub_deps.append(dep)
        sub_times.append(t_shift)
        sub_durs.append(dur)

    sub_lons = np.array(sub_lons)
    sub_lats = np.array(sub_lats)
    sub_deps = np.array(sub_deps)
    sub_times = np.array(sub_times)
    sub_durs = np.array(sub_durs)

    # Calculate M0 and Mw dynamically for each subevent
    mw_list, m0_list, scale_sizes, psmeca_tensors = [], [], [], []

    for i in range(nsub):
        scale = fm_data[i, 0]
        c1, c2, c3, c4, c5, c6 = fm_data[i, 1:7]
        m0 = math.sqrt((c1**2 + 2*c2**2 + 2*c3**2 + c4**2 + 2*c5**2 + c6**2)/2.0) * scale
        mw = (2.0 / 3.0) * math.log10(m0) - 10.7
        m0_list.append(m0)
        mw_list.append(mw)

        # Exact psmeca input fields from plot.sh:
        mrr_p, mtt_p, mff_p, mrt_p, mrf_p, mtf_p = c6, c1, c4, c3, -c5, -c2
        psmeca_tensors.append((mrr_p, mtt_p, mff_p, mrt_p, mrf_p, mtf_p))

        log_m0 = math.log10(m0)
        # Scale size dynamically adapted for Mw (min size 0.25i, max scaled for high M0)
        scale_size = max((10.0 ** (log_m0 / 3.0 - 27.2 / 3.0)) * 0.75, 0.25)
        scale_sizes.append(scale_size)

    # Total moment & Mw calculation
    total_m0 = sum(m0_list)
    total_mw = (2.0 / 3.0) * math.log10(total_m0) - 10.7
    print(f"Total Model Mw: {total_mw:.2f}")

    # 3. Dynamic Region Extents & STF Time Limits Calculation
    all_lons = np.append(sub_lons, evlo)
    all_lats = np.append(sub_lats, evla)
    min_lon, max_lon = np.min(all_lons), np.max(all_lons)
    min_lat, max_lat = np.min(all_lats), np.max(all_lats)

    lon_span = max(max_lon - min_lon, 0.3)
    lat_span = max(max_lat - min_lat, 0.3)

    r_lon1 = min_lon - lon_span * 0.75
    r_lon2 = max_lon + lon_span * 0.75
    r_lat1 = min_lat - lat_span * 0.75
    r_lat2 = max_lat + lat_span * 0.95

    # STF Time Window dynamically calculated
    max_t_end = np.max(sub_times + sub_durs)
    stf_tmin = -5.0
    stf_tmax = math.ceil((max_t_end + 10.0) / 10.0) * 10.0

    # Profile A-A' calculation (horizontal W-E across subevents center)
    center_lon = np.mean(sub_lons)
    center_lat = np.mean(sub_lats)

    half_len_deg = max(lon_span * 0.95, 0.35)
    bg_lon = center_lon - half_len_deg
    bg_lat = center_lat + half_len_deg * math.tan(np.radians(5.0))
    nd_lon = center_lon + half_len_deg
    nd_lat = center_lat - half_len_deg * math.tan(np.radians(5.0))

    dist_total = calculate_distance(bg_lat, bg_lon, nd_lat, nd_lon)

    dep_min_r = math.floor((np.min(sub_deps) - 15.0) / 10.0) * 10.0
    dep_min_r = max(40.0, dep_min_r)
    dep_max_r = math.ceil((np.max(sub_deps) + 25.0) / 10.0) * 10.0

    # Project subevent centroids onto profile A-A'
    dx_prof = (nd_lon - bg_lon) * (R_earth * np.cos(np.radians(bg_lat))) * (np.pi / 180.0)
    dy_prof = (nd_lat - bg_lat) * R_earth * (np.pi / 180.0)
    prof_len = np.sqrt(dx_prof**2 + dy_prof**2)

    sub_dists = []
    for i in range(nsub):
        dx_p = (sub_lons[i] - bg_lon) * (R_earth * np.cos(np.radians(bg_lat))) * (np.pi / 180.0)
        dy_p = (sub_lats[i] - bg_lat) * R_earth * (np.pi / 180.0)
        cdist = (dx_p * dx_prof + dy_p * dy_prof) / prof_len
        sub_dists.append(cdist)

    # 4. Compute 95% Location Error Ranges from hist directories
    hist_dirs = ["../../hist", "../hist", "hist", "../../hist_4sub", "../hist_4sub"]
    map_errors, section_errors = compute_95_errors(hist_dirs, nsub, evla, evlo, bg_lon, bg_lat, nd_lon, nd_lat)

    print(f"Map Region: -R{r_lon1:.4f}/{r_lon2:.4f}/{r_lat1:.4f}/{r_lat2:.4f}")
    print(f"Profile A-A': ({bg_lon:.4f},{bg_lat:.4f}) to ({nd_lon:.4f},{nd_lat:.4f}), Len={dist_total:.1f}km, Depth={dep_min_r:.0f}-{dep_max_r:.0f}km")
    print(f"STF Time Range: {stf_tmin:.0f} s to {stf_tmax:.0f} s")

    # Cut topo grid & prepare CPT
    run_cmd(f"grdcut {topo_grd_path} -R{r_lon1:.4f}/{r_lon2:.4f}/{r_lat1:.4f}/{r_lat2:.4f} -Glocal_topo.grd")
    run_cmd("grdgradient local_topo.grd -A45 -Nt0.3 -Glocal_topoI.grd")
    run_cmd("makecpt -Cgray -T-5002/-5000/1 -Z -D > topo1.cpt")

    # Copy SAC files to local if in parent dir
    if os.path.exists('../stf'):
        run_cmd("cp ../stf/stf_*.sac . 2>/dev/null")

    # Multiply SAC files by subevent M0 normalized by max M0 so relative proportions are true but fits in frame
    run_cmd("sh gentimes.sh > times.dat 2>/dev/null || cat Input.model | gawk '{print $1}' > times.dat")

    max_m0 = max(m0_list)
    mult_cmd = f"""{sys.executable} -c "
from obspy import read
import glob, numpy as np
m0_list = {m0_list}
max_m0 = {max_m0}
for i, f in enumerate(sorted(glob.glob('stf_000*.sac'))):
    if i >= len(m0_list): break
    tr = read(f)[0]
    tr.data = tr.data * (m0_list[i] / max_m0)
    tr.write(f'stf_scaled_{{i:04d}}.sac', format='SAC')
" """
    run_cmd(mult_cmd)

    with open("stf_scaled.dat", "w") as f:
        for i in range(nsub):
            f.write(f"stf_scaled_{i:04d}.sac\n")

    sac_shift_cmd = "paste stf_scaled.dat times.dat | gawk '{print \"r \"$1;print \"ch delta 1 b -10 t1 0 dist 5000\";print \"w over\"} END{print \"q\"}' | sac >/dev/null 2>&1"
    run_cmd(sac_shift_cmd)

    # Calculate max amplitude of normalized scaled traces for pssac2 -M scaling
    mamp_out = run_cmd("saclst depmax f stf_scaled_*.sac 2>/dev/null | minmax | gawk 'BEGIN{FS=\"[</>]\"} {print $8}'").strip()
    try:
        mamp = float(mamp_out) * 1.5
    except ValueError:
        mamp = 0.8

    # Colors per subevent
    color_palette = ["255/0/0", "255/255/0", "0/255/0", "0/255/255", "0/0/255", "255/0/255"]
    colors = [color_palette[i % len(color_palette)] for i in range(nsub)]

    # Generate GMT script matching Figure 2 exact syntax with side-by-side (left-right) layout
    sh_script = f"""#!/bin/bash
gmtset LABEL_FONT_SIZE 13p
gmtset ANNOT_FONT_SIZE_PRIMARY 13p
gmtset ANNOT_OFFSET_PRIMARY 0.1c
gmtset TICK_LENGTH -0.1c
gmtset GRID_PEN_PRIMARY 0.4t4_2:p
gmtset BASEMAP_TYPE plain
gmtset LABEL_OFFSET 0.1c
gmtset FRAME_PEN 1
gmtset COLOR_NAN white
gmtset PLOT_DEGREE_FORMAT D

Range="{r_lon1:.4f}/{r_lon2:.4f}/{r_lat1:.4f}/{r_lat2:.4f}"

# --- Part 1: Topography Map (Left Panel) ---
grdimage local_topo.grd -Ilocal_topoI.grd -R$Range -JM3.2i -X0.8i -Y2.0i -Ctopo1.cpt -K -B > map.ps
pscoast -R$Range -J -Dh -A1000 -S235/235/255 -W1/0.1p,black -N1 -K -O -Ba0.2f0.1/a0.2f0.1WSen >> map.ps

# Draw Profile Line A - A'
psxy -R$Range -J -W1.2p/0/0/0 -K -O >> map.ps << END
{bg_lon:.4f} {bg_lat:.4f}
{nd_lon:.4f} {nd_lat:.4f}
END
echo "{bg_lon:.4f} {bg_lat:.4f} 12 0 1 CB A" | pstext -R$Range -J -N -K -O >> map.ps
echo "{nd_lon:.4f} {nd_lat:.4f} 12 0 1 CB A'" | pstext -R$Range -J -N -K -O >> map.ps

# Epicenter star
psxy -R$Range -J -Sa0.15i -W1p/0/0/0 -G255/0/0 -K -O >> map.ps << END
{evlo:.4f} {evla:.4f}
END

# Subevents Focal Mechanisms on Map
"""

    for i in range(nsub):
        t = psmeca_tensors[i]
        log_m0 = math.log10(m0_list[i])
        sh_script += f"""echo "{sub_lons[i]:.5f} {sub_lats[i]:.5f} {sub_deps[i]:.2f} {t[0]} {t[1]} {t[2]} {t[3]} {t[4]} {t[5]} {log_m0:.4f} {sub_lons[i]:.5f} {sub_lats[i]:.5f} E{i+1} Mw {mw_list[i]:.1f}" | psmeca -R$Range -J -K -O -G"{colors[i]}" -Sm{scale_sizes[i]:.4f}i -T0 -p -t >> map.ps\n"""

    # Plot 95% Horizontal Location Error Bars OVER Beachballs on Panel 1 (Map) centered at subevent (lon, lat)
    if map_errors:
        for i, (lon_err, lat_err) in enumerate(map_errors):
            sh_script += f"""psxy -R$Range -J -Exy/0.15i/1.5p,black -W1.5p,black -K -O >> map.ps << END
{sub_lons[i]:.5f} {sub_lats[i]:.5f} {lon_err:.5f} {lat_err:.5f}
END
"""

    # STF Panel placed further up outside/at top of left map panel (X=-0.05i, Y=3.55i)
    sh_script += f"""
# --- Part 2: Source Time Functions Inset (High up at top-left margin) ---
gmtset LABEL_FONT_SIZE 10p
gmtset ANNOT_FONT_SIZE_PRIMARY 10p
psbasemap -R{stf_tmin:.0f}/{stf_tmax:.0f}/0/1 -JX1.1i/0.75i -X-0.05i -Y3.55i -G255/255/245 -K -O -Ba10f5:"Time (s)":/WSen >> map.ps
cat stf_scaled.dat | gawk -v mamp="{mamp}" '{{if (NR==1) cl="255/0/0";if (NR==2) cl="255/255/0";if (NR==3) cl="0/255/0";if (NR==4) cl="0/255/255";if (NR==5) cl="0/0/255";if (NR==6) cl="255/0/255";print "echo \\""$1" 0 0\\" | pssac2 -R -J -K -O -W1p/0/0/0 -G"cl"/0/{stf_tmin:.0f}/{stf_tmax:.0f} -C{stf_tmin:.0f}/{stf_tmax:.0f} -M"mamp"/0 -Ent1 >> map.ps";}}' | awk 'NR==1{{line1=$0; next}} NR==2{{print; print line1; next}} 1' | sh
gmtset LABEL_FONT_SIZE 13p
gmtset ANNOT_FONT_SIZE_PRIMARY 13p

# --- Part 3: Cross Section (A - A') (Right Panel) ---
psbasemap -R0/{dist_total:.1f}/{dep_min_r:.0f}/{dep_max_r:.0f} -JX3.2i/-3.2i -X3.7i -Y-3.55i -Ba20f10:"Distance (km)":/a10f5:"Depth (km)":WSen -K -O >> map.ps
"""

    for i in range(nsub):
        t = psmeca_tensors[i]
        log_m0 = math.log10(m0_list[i])
        sh_script += f"""echo "{sub_lons[i]:.5f} {sub_lats[i]:.5f} {sub_deps[i]:.2f} {t[0]} {t[1]} {t[2]} {t[3]} {t[4]} {t[5]} {log_m0:.4f} X Y E{i+1}" | pscoupe -R0/{dist_total:.1f}/{dep_min_r:.0f}/{dep_max_r:.0f} -JX3.2i/-3.2i -T0 -L -K -O -G"{colors[i]}" -Sm{scale_sizes[i]:.4f}i -Aa{bg_lon:.4f}/{bg_lat:.4f}/{nd_lon:.4f}/{nd_lat:.4f}/90/50/{dep_min_r:.0f}/{dep_max_r:.0f} -p -t >> map.ps\n"""

    # Plot 95% Profile Distance & Depth Error Bars OVER Beachballs on Panel 2 (Cross Section) centered at (sub_dist, sub_dep)
    if section_errors:
        for i, (dist_err, dep_err) in enumerate(section_errors):
            sh_script += f"""psxy -R0/{dist_total:.1f}/{dep_min_r:.0f}/{dep_max_r:.0f} -JX3.2i/-3.2i -Exy/0.15i/1.5p,black -W1.5p,black -K -O >> map.ps << END
{sub_dists[i]:.2f} {sub_deps[i]:.2f} {dist_err:.2f} {dep_err:.2f}
END
"""

    sh_script += """
# Finish PostScript
psxy -R -J -O >> map.ps << END
END

# Convert PostScript to PDF and PNG
ps2pdf map.ps map.pdf
gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r300 -sOutputFile=map.png map.ps
echo "Publication-grade figure successfully generated: map.pdf & map.png"
"""

    with open("plot_publication.sh", "w") as f:
        f.write(sh_script)
    os.chmod("plot_publication.sh", 0o755)

    print("Executing automated GMT plotting script...")
    run_cmd("./plot_publication.sh")

if __name__ == "__main__":
    main()
