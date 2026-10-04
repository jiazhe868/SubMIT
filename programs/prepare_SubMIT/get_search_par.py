import sys

def generate_search_par_file(event_name, mw, nsubeve):
    # Define duration; giant events (M >~ 8) use length-based duration
    # (same rule as gen_par.py)
    duration = round(5 * 10**(0.5 * (mw - 6)))
    dur_len = round(10 ** (0.59 * mw - 2.44) / 2.5)
    if dur_len > 1.5 * duration:
        print(f"get_search_par: length-based duration {dur_len} s (giant event)")
        duration = dur_len

    # Set burn_in and nsample
    burn_in = nsubeve * 1000
    nsample = nsubeve * 1000

    # Set neq_min and neq_max
    neq_min = nsubeve
    neq_max = nsubeve

    # Read x_min, x_max, y_min, y_max from edges.dat (nonzero aftershock-density
    # region written by caldens.py)
    with open('edges.dat', 'r') as edges_file:
        edges = edges_file.readline().split()
        x_min = float(edges[0])
        x_max = float(edges[1])
        y_min = float(edges[2])
        y_max = float(edges[3])
    # always allow at least +/-15 km around the hypocenter, but never beyond the
    # density grid: finv's seisprior() returns 0 outside seisdens.dat, so any
    # wider search range would be silently ineffective
    try:
        with open('seisgrids.dat') as g:
            gvals = [float(v) for v in g.read().split()]
        gx_min, gx_max, gy_min, gy_max = gvals[0], gvals[1], gvals[2], gvals[3]
    except Exception:
        gx_min, gx_max, gy_min, gy_max = x_min, x_max, y_min, y_max
    x_min = max(min(x_min, -15.0), gx_min)
    x_max = min(max(x_max, 15.0), gx_max)
    y_min = max(min(y_min, -15.0), gy_min)
    y_max = min(max(y_max, 15.0), gy_max)

    # Read z_min, z_max from distdep.dat
    with open('distdep.dat', 'r') as distdep_file:
        distdep = distdep_file.readline().split()
        z_min = float(distdep[3])
        z_max = float(distdep[4])
    # CRUSTAL events stay in the crust: for a continental hypocenter above the
    # Moho, cap the depth search at Moho+4 km (aftershock catalogs often carry
    # fixed 10 km depths, so the prior cannot do this job). Intraslab events
    # (hypocenter below Moho, e.g. Calama 109 km) keep the full evdp+/-20 range.
    try:
        import glob as _glob
        vmf = (_glob.glob(f"../gf_{event_name}/loc/vmodel.txt")
               or _glob.glob("../gf_*/loc/vmodel.txt"))[0]
        moho = sum(float(l.split()[3]) for l in open(vmf)
                   if l.split() and float(l.split()[0]) < 7.5)
        evdp = float(open('mainshock.dat').read().split()[2])
        # only for STRIKE-SLIP events (B-plunge >= 60 from the trusted total
        # tensor below): a megathrust INTERFACE event ruptures below the local
        # Moho (Kamchatka M8.8: depths to 36 km vs offshore Moho ~25)
        v = [float(x) for x in open('totalmt.dat').read().split()[:5]]
        import numpy as _np2
        _M = _np2.array([[v[0], v[2], v[3]], [v[2], v[1], v[4]],
                         [v[3], v[4], -(v[0] + v[1])]])
        _w, _vec = _np2.linalg.eigh(_M)
        _b = _vec[:, _np2.argsort(_w)[1]]
        _bplunge = _np2.degrees(_np2.arcsin(abs(_b[2])))
        src_ok = open('totalmt_source.txt').read().strip() in (
            "gcmt-local", "usgs", "usgs-sum")
        if src_ok and _bplunge >= 60 and 10.0 < moho and evdp <= moho and z_max > moho + 4:
            print(f"get_search_par: crustal event (evdp {evdp:.0f} km <= Moho "
                  f"{moho:.0f} km) -> z_max capped {z_max:.0f} -> {moho + 4:.0f}")
            z_max = round(moho + 4)
    except Exception:
        pass

    # Set dura_min and dura_max; never let the magnitude scaling squeeze the
    # subevent duration below its own estimate of the source duration
    dura_min = 3
    if nsubeve == 1:
        dura_max = 30
    elif nsubeve == 2:
        dura_max = 20
    else:
        dura_max = 15
    # cap at 35 s: a "subevent" is one coherent pulse - the magnitude scaling
    # (M8.8 -> 131 s) is the WHOLE source duration, not a subevent's
    # (Kamchatka M8.8 production used dura_max=30)
    dura_max = min(35, max(dura_max, round(duration) + 5))

    # Set cen_min and cen_max. cen_max = duration + 6: benchmarked 2026-07-29 -
    # finv staggers the INITIAL subevent centroids across [cen_min, cen_max], so
    # an over-generous cen_max (e.g. 1.5*duration) forces late subevents to start
    # in the coda window and traps them there (broke misfit-vs-nsub monotonicity;
    # 3sub best 0.5686 at cen_max 26 vs 0.4889 at cen_max 20 on the Calama M6.9).
    cen_min = 1.5
    cen_max = round(duration) + 6
    # extended-source / doublet widening (same test as gen_par.py): if the
    # aftershock footprint reaches far beyond the mainshock's rupture length,
    # late moment release (e.g. a doublet partner) must be reachable in time -
    # r_corner at a conservative 2.5 km/s apparent speed, plus the duration
    # (2026 Venezuela M7.2: partner M7.5 fit at cen ~99 s, 172 km east)
    L_rup = 10 ** (0.59 * mw - 2.44)     # Wells & Coppersmith (1994) SRL, km
    try:
        e = [float(v) for v in open('edges.dat').read().split()[:4]]
        r_corner = (max(abs(e[0]), abs(e[1])) ** 2
                    + max(abs(e[2]), abs(e[3])) ** 2) ** 0.5
    except Exception:
        r_corner = 0.0
    # kinematic reachability: the far corner of the searched region must be
    # reachable in time at ~2.5 km/s - covers BOTH doublets (Venezuela: partner
    # at 172 km / 99 s) and giant single ruptures whose true duration exceeds
    # the magnitude scaling (Kamchatka M8.8: 500 km / 154 s vs scaling 126 s).
    # Trigger only when traversal time clearly exceeds the scaled duration so
    # ordinary events (Calama) keep the validated duration+6.
    if r_corner > max(1.5 * L_rup, L_rup + 40.0):
        cen_max = int(max(cen_max, round(r_corner / 2.5 + dura_max + 6)))
        # EXTENDED/doublet sources (user physics 2026-08-12): cascading
        # high-frequency asperity breaks (short subevents, fit by P) can ride
        # on a SMOOTH slow rupture pronounced in long-period SH (band reaches
        # 0.005 Hz). The standard dura_max cap makes that smooth component
        # unrepresentable (2026 Venezuela: durations railed at the cap, SH
        # underfit) - give one very-long-duration subevent headroom.
        dura_ext = int(min(90, max(dura_max, round(0.5 * cen_max))))
        if dura_ext > dura_max:
            print(f"get_search_par: extended source - dura_max {dura_max} "
                  f"-> {dura_ext} (smooth-rupture component headroom)")
            dura_max = dura_ext
        print(f"get_search_par: footprint {r_corner:.0f} km (L {L_rup:.0f} km, "
              f"traversal {r_corner/2.5:.0f} s vs duration {duration:.0f} s) "
              f"-> cen_max {cen_max}")

    # Other fixed values
    vr_min = 0.5
    vr_max = 3.5
    theta_min = 0
    theta_max = 360

    # Create search_par.file with the specified format
    with open('search_par.file', 'w') as par_file:
        par_file.write(f"burn_in {burn_in}\n")
        par_file.write(f"nsample {nsample}\n")
        par_file.write(f"neq_min {neq_min}\n")
        par_file.write(f"neq_max {neq_max}\n")
        par_file.write(f"x_min {x_min}\n")
        par_file.write(f"x_max {x_max}\n")
        par_file.write(f"y_min {y_min}\n")
        par_file.write(f"y_max {y_max}\n")
        par_file.write(f"z_min {z_min}\n")
        par_file.write(f"z_max {z_max}\n")
        par_file.write(f"dura_min {dura_min}\n")
        par_file.write(f"dura_max {dura_max}\n")
        par_file.write(f"vr_min {vr_min}\n")
        par_file.write(f"vr_max {vr_max}\n")
        par_file.write(f"theta_min {theta_min}\n")
        par_file.write(f"theta_max {theta_max}\n")
        par_file.write(f"cen_min {cen_min}\n")
        par_file.write(f"cen_max {cen_max}\n")
        # Rupture-causality cap: negative = factor of Vs inside finv.
        # Default 1.0*Vs; shallow strike-slip events may rupture supershear, so
        # allow up to 1.5*Vs when depth <= 35 km AND the GCMT total mechanism is
        # strike-slip (null/B axis plunge >= 60 deg). Falls back to 1.0*Vs when
        # totalmt.dat/mainshock.dat are unavailable.
        vfac = -1.0
        try:
            import numpy as _np
            # only trust the mechanism if totalmt.dat is a REAL solution for
            # this event (gcmt-local/usgs); a stale template must never grant
            # supershear headroom
            src = open('totalmt_source.txt').read().strip()
            assert src in ("gcmt-local", "usgs", "usgs-sum"), f"totalmt source: {src}"
            dep = float(open('mainshock.dat').read().split()[2])
            v = [float(x) for x in open('totalmt.dat').read().split()[:5]]
            mxx, myy, mxy, mxz, myz = v
            M = _np.array([[mxx, mxy, mxz], [mxy, myy, myz],
                           [mxz, myz, -(mxx + myy)]])
            w, vec = _np.linalg.eigh(M)
            b_axis = vec[:, _np.argsort(w)[1]]          # intermediate (null) axis
            plunge = _np.degrees(_np.arcsin(abs(b_axis[2])))
            if dep <= 35 and plunge >= 60:
                vfac = -1.5
                print(f"get_search_par: shallow strike-slip (dep {dep:.0f} km, "
                      f"B-plunge {plunge:.0f} deg) -> vmax_cau 1.5*Vs")
        except Exception:
            pass
        par_file.write(f"vmax_cau {vfac}\n")  # negative = factor of Vs (finv)

if __name__ == "__main__":
    event_name = sys.argv[1]
    mw = float(sys.argv[2])
    nsubeve = int(sys.argv[3])
    generate_search_par_file(event_name, mw, nsubeve)

