import sys
import numpy as np
from scipy.interpolate import interp1d

def read_weights(file_path):
    """Read weights from the weights.dat file."""
    with open(file_path, 'r') as f:
        weights = f.readline().strip().split()
        return float(weights[0]), float(weights[1]), float(weights[2])

def read_distdep(file_path):
    """Read distance and depth parameters from the distdep.dat file."""
    with open(file_path, 'r') as f:
        values = f.readline().strip().split()
        return int(values[0]), int(values[1]), int(values[2]), int(values[3]), int(values[4]), int(values[5])

def count_rows(file_path):
    """Count the number of rows in a file."""
    with open(file_path, 'r') as f:
        return sum(1 for _ in f)

def load_fine_iasp91():
    """
    Loads a finer IASP91 model data for depth, Vp, and Vs.
    """
    # Finer depth points for IASP91 model (depth in km, Vp and Vs in km/s)
    depth = np.array([0, 10, 20, 35, 50, 70, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700])
    vp = np.array([5.80, 6.20, 6.50, 8.04, 8.10, 8.15, 8.04, 8.17, 8.29, 8.35, 8.50, 8.75, 9.03, 9.10, 9.36, 9.48, 9.60, 9.60, 9.60])
    vs = np.array([3.36, 3.60, 3.85, 4.47, 4.55, 4.62, 4.47, 4.63, 4.74, 4.80, 4.87, 4.92, 5.08, 5.15, 5.30, 5.42, 5.52, 5.52, 5.52])
    return depth, vp, vs

def calculate_velocities_at_depth(depth):
    """
    Given a source depth, calculates Vp and Vs using a finer 1D velocity model (IASP91).
    """
    # Load the finer velocity model data
    depth_model, vp_model, vs_model = load_fine_iasp91()
    
    # Create interpolation functions
    vp_interp = interp1d(depth_model, vp_model, kind='linear', fill_value="extrapolate")
    vs_interp = interp1d(depth_model, vs_model, kind='linear', fill_value="extrapolate")
    
    # Calculate interpolated Vp and Vs at the given depth
    vp = vp_interp(depth)
    vs = vs_interp(depth)
    vp = float(f"{vp:.1f}")
    vs = float(f"{vs:.1f}")
    return vp, vs

def read_deltat():
    """Central dt from programs/submit.conf (default 1 s)."""
    import os
    for c in ("../../programs/submit.conf", "../programs/submit.conf"):
        if os.path.exists(c):
            for line in open(c):
                if line.strip().startswith("DELTAT="):
                    return float(line.strip().split("=")[1])
    return 1.0

def generate_par_file(event_name, mw, evlo, evla, evdp):
    """Generate Par.file with specified parameters and formatting."""
    # Calculate duration; for GIANT events (M >~ 8) the amplitude scaling
    # underestimates the true source duration - use the Wells-Coppersmith
    # rupture length traversed at 2.5 km/s instead (Kamchatka M8.8: scaling
    # 126 s vs actual ~200 s rupture; production windows were sized for 280 s)
    duration = 5 * 10 ** (0.5 * (mw - 6))
    dur_len = 10 ** (0.59 * mw - 2.44) / 2.5
    if dur_len > 1.5 * duration:
        print(f"gen_par: length-based duration {dur_len:.0f} s replaces "
              f"magnitude scaling {duration:.0f} s (giant event)")
        duration = dur_len
    
    # extended-source / doublet detection: aftershock footprint (edges.dat,
    # written by caldens.py) far beyond the mainshock's own rupture length
    # means moment release the duration-scaled windows would truncate
    # (2026 Venezuela: M7.5 partner 30 s later / 140 km east of the M7.2)
    L_rup = 10 ** (0.59 * mw - 2.44)     # Wells & Coppersmith (1994) SRL, km
    r_corner = 0.0
    try:
        e = [float(v) for v in open('edges.dat').read().split()[:4]]
        r_corner = (max(abs(e[0]), abs(e[1])) ** 2
                    + max(abs(e[2]), abs(e[3])) ** 2) ** 0.5
    except Exception:
        pass
    extended = r_corner > max(1.5 * L_rup, L_rup + 40.0)
    # persist the decision: step3's tighten_bounds_from_1sub.sh must NOT
    # shrink cen_max for extended/doublet sources (the 1-sub scan only
    # captures the first event of a doublet)
    with open('extended.flag', 'w') as ef:
        ef.write('1' if extended else '0')
    # widen so the latest reachable subevent (same kinematic cen_max estimate
    # as get_search_par: traversal of the search region at 2.5 km/s) plus its
    # duration fits in the body windows. Fires for doublets (Venezuela: 152 s
    # ~ hand-tuned 150) AND giant single ruptures whose true duration exceeds
    # the magnitude scaling (Kamchatka M8.8: ~280 s vs scaling 126 s)
    kinematic = False   # giant-rupture case handled by length-based duration
    dura_est = min(35, round(duration) + 5)
    cen_est = round(r_corner / 2.5 + dura_est + 6)
    t_extra = 0
    if extended or kinematic:
        t_extra = max(0, cen_est + dura_est + 15 - round(duration + 70))
        print(f"gen_par: footprint {r_corner:.0f} km (L {L_rup:.0f} km, "
              f"traversal {r_corner/2.5:.0f} s vs duration {duration:.0f} s) "
              f"-> body windows widened by {t_extra} s")
        if extended:
            with open("extended_source.txt", "w") as f:
                f.write(f"{r_corner:.0f} {L_rup:.0f}\n")

    # frequency bands: the body-wave corner does NOT scale with magnitude
    # (M7.7 Myanmar production uses 0.2 Hz) - it drops to 0.1 Hz only for a
    # genuine multi-mainshock compound source, which gen_totalmt detects from
    # the catalog and patches into Par.file (2026 Venezuela doublet: 0.1 Hz).
    hf_body = 0.2
    # Rayleigh corner stays near 0.1 Hz even for large events - lowering it
    # discards temporal resolution (calibration: M6.9 Calama 0.15, M7.7
    # Myanmar 0.13; the 0.05 once used on Venezuela was too low)


    # Set parameters for time windows
    bg_timeP = -20
    bg_timeSH = -20
    bg_timeRayl = -20
    # intermediate-depth events (user rule 2026-08-13): depth phases arrive
    # ~2z/v after the direct phase (sS-S ~ 2*120/4 = 60 s at 120 km) and are
    # the depth-resolving signal - windows must contain them. The +70 already
    # holds ~30 s of shallow phase allowance; add the excess for deep sources.
    t_depth = round(max(0.0, 2.0 * evdp / 4.0 - 30.0)) if evdp > 60 else 0
    nd_timeP = round(duration + 70) + t_extra + t_depth
    nd_timeSH = round(duration + 70) + t_extra + t_depth
    # STF window must contain cen_min - dura_max/2 (subevent onset) at the
    # low end and the latest reachable subevent at the high end
    stf_btime = min(-10, -(dura_est // 2) - 2)
    # STF window must contain the latest reachable subevent entirely
    stf_etime = (cen_est + dura_est + 5) if (extended or kinematic) else round(duration + 30)
    deltat = read_deltat()
    vp, vs = calculate_velocities_at_depth(evdp)
    # Count the number of rows (stations) in each file
    num_sta_P = count_rows('stations.info')
    num_sta_SH = count_rows('stationsSH.info')
    num_sta_Rayl = count_rows('stationsloc.info')
    
    # Read weights from weights.dat
    weightP, weightSH, weightRayl = read_weights('weights.dat')
    
    # Read distance and depth parameters from distdep.dat
    bg_localdist, nd_localdist, localdist_interval, bg_depth, nd_depth, depth_interval = read_distdep('distdep.dat')

    # farthest ACTUAL station (nd_localdist carries a +200 km GF-grid buffer
    # and a 950 km cap, so it over/under-states the real aperture)
    import math
    max_sta = 0.0
    try:
        for line in open('stationsloc.info'):
            t = line.split()
            if len(t) >= 3:
                dx = (float(t[1]) - evlo) * 111.32 * math.cos(math.radians(evla))
                dy = (float(t[2]) - evla) * 110.574
                max_sta = max(max_sta, math.hypot(dx, dy))
    except Exception:
        pass
    if max_sta <= 0:
        max_sta = max(nd_localdist - 200, 0)

    # Rayleigh corner is APERTURE-limited, not magnitude-limited: 1D-model
    # coherence dies with path length. Calibrated on production events:
    # Calama 390 km -> 0.15 works; Myanmar 570 km -> 0.13 works;
    # Venezuela 731-885 km -> 0.106 destroys the joint MT solve, 0.05 works.
    # Piecewise through (450,0.15) (600,0.13) (750,0.06), floor 0.05.
    if max_sta <= 450:
        hf_rayl = 0.15
    elif max_sta <= 570:
        hf_rayl = round(0.15 - 0.02 * (max_sta - 450) / 120.0, 3)
    else:
        # steep: 0.069 already halves the MT amplitudes on 731-km paths
        hf_rayl = round(max(0.05, 0.13 - 0.08 * (max_sta - 570) / 130.0), 3)

    # Rayleigh window scales with the true aperture (a fixed 210 s cut off
    # surface waves beyond ~550 km); 3.0 km/s group speed + 30 s margin keeps
    # the validated 210 s for <=450 km networks; clamped to maxnpts=500
    nd_timeRayl = max(210, round(max_sta / 3.0 + duration + 30)) + t_extra
    # maxnpts=500 in sub_header.h; npts spans both endpoints, hence the -1
    nd_max = bg_timeRayl + int((500 - 1) * deltat)
    if nd_timeRayl > nd_max:
        print(f"gen_par: nd_timeRayl {nd_timeRayl} clamped to {nd_max} "
              f"(maxnpts=500; enlarge sub_header.h to extend)")
        nd_timeRayl = nd_max

    # Local max time shift scales with the event's true aperture: ~10%
    # velocity heterogeneity at ~3 km/s over the farthest local path
    # (450 km -> 15 s, 100 km -> 3 s), floored at 3 s
    maxshft_rayl = int(min(20, max(3, round(0.10 * max_sta / 3.0))))

    # Set green function character format
    greenchar1 = f"../gf_{event_name}/greenFuncDir_disp/vmodel_"
    greenchar2 = "/"
    greenchar3 = ".grn."
    
    # Generate Par.file content
    par_content = f"""\
bg_timeP= {bg_timeP} #begining time of inverse time window (P)
nd_timeP= {nd_timeP} #ending time of inverse time window (P)
bg_timeSH= {bg_timeSH} #begining time of inverse time window (SH)
nd_timeSH= {nd_timeSH} #ending time of inverse time window (SH)
bg_timeRayl= {bg_timeRayl} #begining time of inverse time window (SH)
nd_timeRayl= {nd_timeRayl} #ending time of inverse time window (SH)
stf_btime= {stf_btime} #begining time of source time function time window
stf_etime= {stf_etime} #ending time of source time function time window
deltat= {deltat} #sampling rate of data (from programs/submit.conf)
num_sta_P= {num_sta_P} #Number of all stations (P)
num_sta_SH= {num_sta_SH} #Number of all stations (SH)
num_sta_Rayl= {num_sta_Rayl} #Number of all stations (P)
stainfofileP= stations.info #Station distance (in degree) and azimuth (P)
stainfofileSH= stationsSH.info #Station distance (in degree) and azimuth (SH)
stainfofileRayl= stationsloc.info
weightP= {weightP} #Weight for P waves in the inversion (0->do not use)
weightSH= {weightSH} #Weight for SH waves in the inversion (0->do not use)
weightRayl= {weightRayl}
InputModelFile= Input.model #Input model parameters
num_subevent= 6 #Number of subevents in input model
Low_frequencyP= 0.005 #Lower frequency bound in Hz (P)
High_frequencyP= {hf_body} #Lower frequency bound in Hz (P)
Low_frequencySH= 0.005 #Lower frequency bound in Hz (SH)
High_frequencySH= {hf_body} #Lower frequency bound in Hz (SH)
Low_frequencyRayl= 0.02 #Lower frequency bound in Hz (P)
High_frequencyRayl= {hf_rayl} #Lower frequency bound in Hz (P)
Tikhonov_alpha= 5e-4 #Parameter for Tikhonov regularization, try using values between 1e-4 to 1e-3 first
CMTscaling= 1e-9 #Parameter for CMT scaling, ~ >1e-5: start to have constraints on total MT
evlo= {evlo}
evla= {evla}
evvp= {vp} #Average Vp near the source depth
evvs= {vs}  #Average Vs near the source depth
DCconstrain= 1 #DC constrain: 2-weak, 0.5-strong
bg_localdist= {bg_localdist}  # begining distance for local stations (km)
nd_localdist= {nd_localdist} # ending distance for local stations (km)
localdist_interval= {localdist_interval}  # distace interval for local stations (km)
bg_depth= {bg_depth} # begining depth in the search (km).
nd_depth= {nd_depth} # ending depth in the search (km).
depth_interval= 2 # depth interval in the seach (km).
maxshftP= 2 # max allowable time shift (s) for P cross-correlation alignment
maxshftSH= 5 # max allowable time shift (s) for SH
maxshftRayl= {maxshft_rayl} # max shift (s) for local waves: 10% velocity error x max distance / 3 km/s
greenchar1= {greenchar1} # "greenchar" series are the format of green's function paths. For example, if the path was "greenFuncDir/vmodel_20/5871.grn." for depth=20km and dist=5871km, then greenchar1="greenFuncDir/vmodel_", greenchar2="/", greenchar3=".grn."
greenchar2= {greenchar2} # See details for greenchar1
greenchar3= {greenchar3} # See details for greenchar1
TotalMomentTensor= totalmt.dat #moment tensor summation should be consistent with point source solution. Unit is 1e27 dyne-cm
"""

    # Write to Par.file
    with open('Par.file', 'w') as f:
        f.write(par_content)

if __name__ == "__main__":
    # Ensure required arguments are passed
    if len(sys.argv) != 6:
        print("Usage: python gen_par.py <event_name> <mw> <evlo> <evla> <evdp>")
        sys.exit(1)

    # Extract arguments
    event_name = sys.argv[1]
    mw = float(sys.argv[2])
    evlo = float(sys.argv[3])
    evla = float(sys.argv[4])
    evdp = float(sys.argv[5])
    
    # Generate Par.file
    generate_par_file(event_name, mw, evlo, evla, evdp)
    print("Par.file generated successfully.")

