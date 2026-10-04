"""Generate totalmt.dat (five rows: Mxx, Myy, Mxy, Mxz, Myz, unit 1e27 dyne-cm)
for the total-moment-tensor constraint, from GCMT when available.

Search order: local catalog programs/mt_SubMIT/1990_2023_gcmt.dat (cols: lon lat
dep mrr mtt mpp mrt mrp mtp exponent X Y ID[MMDDYY?]), else USGS online
moment-tensor product. RTP -> XYZ per mrtp2xyz.sh (x north, y east, z down):
Mxx=Mtt, Myy=Mpp, Mxy=-Mtp, Mxz=Mrt, Myz=-Mrp.
On any failure the existing totalmt.dat template is kept (harmless while
CMTscaling ~ 0 in Par.file) and a loud warning is printed.

usage: python gen_totalmt.py <isotime> <evlo> <evla> <evdp> <mw>
"""
import os
import sys
from datetime import datetime

def mark(source):
    with open("totalmt_source.txt", "w") as f:
        f.write(source + "\n")
    # any TRUSTED total MT enables the total-moment-tensor constraint: the
    # band-limited misfit does not anchor the absolute moment, and without
    # the constraint the total magnitude drifts (Kamchatka M8.8 -> 8.34,
    # Calama M6.9 -> 7.12). The constraint rows are SCALE-INVARIANT in
    # sub_forward.c (weighted by data-norm/target-moment), so one
    # dimensionless strength works for every event. CMTscaling ~ 1/tolerance:
    # lam=3 tolerates ~33% moment deviation (dMw ~ 0.1) before the constraint
    # costs as much as the data misfit - a soft anchor that recovers the
    # target in SAMPLING (kinematics adapt; Calama sampled identically at 3
    # and 5) while still letting a genuinely different total be reported.
    # A template tensor (wrong event) keeps the constraint OFF (1e-9).
    import os
    import re
    if source in ("gcmt-local", "usgs", "usgs-sum") and os.path.exists("Par.file"):
        txt = open("Par.file").read()
        txt2 = re.sub(r"CMTscaling= \S+", "CMTscaling= 3", txt)
        if txt2 != txt:
            open("Par.file", "w").write(txt2)
            print(f"gen_totalmt: trusted tensor ({source}) -> CMTscaling= 3 (scale-invariant)")

def write_totalmt(mtt, mpp, mtp, mrt, mrp, scale):
    rows = [mtt * scale, mpp * scale, -mtp * scale, mrt * scale, -mrp * scale]
    with open("totalmt.dat", "w") as f:
        f.write("\n".join(f"{v:.4f}" for v in rows) + "\n")
    print(f"gen_totalmt: wrote totalmt.dat (Mxx Myy Mxy Mxz Myz = "
          + " ".join(f"{v:.3f}" for v in rows) + " x1e27 dyne-cm)")

def from_local(t0, lon, lat):
    for c in ("../../programs/mt_SubMIT/1990_2023_gcmt.dat",
              "../programs/mt_SubMIT/1990_2023_gcmt.dat"):
        if not os.path.exists(c):
            continue
        best = None
        for line in open(c):
            p = line.split()
            if len(p) < 13:
                continue
            evid = p[12]
            try:
                d = datetime.strptime(evid[:6], "%m%d%y")
            except ValueError:
                continue
            if d.date() != t0.date():
                continue
            dist = ((float(p[0]) - lon) ** 2 + (float(p[1]) - lat) ** 2) ** 0.5
            if dist < 1.5 and (best is None or dist < best[0]):
                best = (dist, p)
        if best:
            p = best[1]
            mrr, mtt, mpp, mrt, mrp, mtp = (float(v) for v in p[3:9])
            scale = 10.0 ** (int(p[9]) - 27)
            write_totalmt(mtt, mpp, mtp, mrt, mrp, scale)
            mark("gcmt-local")
            return True
    return False

def fetch_tensor(detail_url):
    import requests
    detail = requests.get(detail_url, timeout=60).json()
    for prod in detail["properties"]["products"].get("moment-tensor", []):
        pr = prod["properties"]
        try:
            return [float(pr[f"tensor-{k}"])
                    for k in ("mrr", "mtt", "mpp", "mrt", "mrp", "mtp")]
        except KeyError:
            continue
    return None

def from_usgs(t0, lon, lat, mw):
    import requests
    from datetime import timedelta, timezone
    # +/-2 min window: catalog origin times can differ by seconds between
    # agencies (the 2026 Venezuela M7.2 is 1.2 s EARLIER in the USGS catalog
    # than in our event_list, which a [t0, t0+2min] window missed)
    tq = t0.replace(microsecond=0)
    # extended sources (marker from gen_par.py): the inverted waveforms contain
    # every large event in the sequence, so the total-MT constraint must be the
    # SUM of their tensors, and the search box must span the whole footprint
    # (the 2026 Venezuela M7.5 partner is 1.3 deg east of the M7.2)
    extended = os.path.exists("extended_source.txt")
    dbox = 1.0
    if extended:
        try:
            dbox += float(open("extended_source.txt").read().split()[0]) / 111.0
        except Exception:
            dbox += 2.0
    q = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
         f"&starttime={(tq - timedelta(minutes=2)).isoformat()}"
         f"&endtime={(tq + timedelta(minutes=2)).isoformat()}"
         f"&minlatitude={lat-dbox}&maxlatitude={lat+dbox}"
         f"&minlongitude={lon-dbox}&maxlongitude={lon+dbox}"
         f"&minmagnitude={mw-0.6}&producttype=moment-tensor")
    r = requests.get(q, timeout=60)
    r.raise_for_status()
    feats = r.json().get("features", [])
    if not feats:
        return False
    # closest match in time+space first (a doublet partner can share the window)
    def score(f):
        p, (flo, fla) = f["properties"], f["geometry"]["coordinates"][:2]
        dt = abs(p["time"] / 1000.0
                 - t0.replace(tzinfo=timezone.utc).timestamp())
        return dt + 100.0 * ((flo - lon) ** 2 + (fla - lat) ** 2) ** 0.5
    feats.sort(key=score)
    use = feats if extended else feats[:1]
    total, names = None, []
    for f in use:
        ten = fetch_tensor(f["properties"]["detail"])
        if ten is None:
            continue
        total = ten if total is None else [a + b for a, b in zip(total, ten)]
        names.append(f["properties"]["title"])
    if total is None:
        return False
    if len(names) > 1:
        print(f"gen_totalmt: extended source - summed {len(names)} tensors: "
              + " + ".join(names))
        # a true multi-mainshock compound source: overlapping ruptures are not
        # fittable at 0.2 Hz with point-source subevents - drop the body-wave
        # corner to 0.1 Hz (single long ruptures like Myanmar M7.7 keep 0.2)
        import re
        if os.path.exists("Par.file"):
            txt = open("Par.file").read()
            txt2 = re.sub(r"High_frequencyP= \S+", "High_frequencyP= 0.1", txt)
            txt2 = re.sub(r"High_frequencySH= \S+", "High_frequencySH= 0.1", txt2)
            if txt2 != txt:
                open("Par.file", "w").write(txt2)
                print("gen_totalmt: doublet -> body-wave corner 0.1 Hz")
    mrr, mtt, mpp, mrt, mrp, mtp = total
    write_totalmt(mtt, mpp, mtp, mrt, mrp, 1e7 / 1e27)  # N-m -> 1e27 dyne-cm
    mark("usgs" if len(names) == 1 else "usgs-sum")
    return True

if __name__ == "__main__":
    # SUBMIT_CATALOG_CACHE=<dir>: replay the saved tensor + its provenance
    cache = os.environ.get("SUBMIT_CATALOG_CACHE")
    if cache and os.path.exists(os.path.join(cache, "totalmt.dat")):
        import shutil
        shutil.copy(os.path.join(cache, "totalmt.dat"), "totalmt.dat")
        src = "gcmt-local"
        if os.path.exists(os.path.join(cache, "totalmt_source.txt")):
            src = open(os.path.join(cache, "totalmt_source.txt")).read().strip()
        mark(src)
        print(f"gen_totalmt: replayed cached tensor ({src}) from {cache}")
        sys.exit(0)
    iso, lon, lat = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    mw = float(sys.argv[5]) if len(sys.argv) > 5 else 6.0
    for fmt in ("%Y-%m-%dT%H%M%S", "%Y/%m/%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S"):
        try:
            t0 = datetime.strptime(iso, fmt)
            break
        except ValueError:
            t0 = None
    try:
        ok = t0 is not None and (from_local(t0, lon, lat) or from_usgs(t0, lon, lat, mw))
    except Exception as exc:
        print(f"gen_totalmt: lookup failed ({exc})")
        ok = False
    if not ok:
        mark("template-UNRELIABLE")
        print("WARNING gen_totalmt: no GCMT/USGS moment tensor found - keeping the "
              "totalmt.dat TEMPLATE (WRONG event!). Do NOT enable CMTscaling "
              "without fixing totalmt.dat manually.")
