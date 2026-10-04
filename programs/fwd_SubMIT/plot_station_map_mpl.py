#!/usr/bin/env python3
"""Station map of an inversion (run in a fwd_<event>_<N>sub directory).

Left: teleseismic stations on an azimuthal-equidistant map centred on the
epicenter (rings every 30 degrees), marked by the wave types they contribute.
Right: regional stations around the epicenter.
Reads Par.file (evlo/evla), stations.info (P), stationsSH.info (SH),
stationsloc.info (regional) and the SAC headers in data/.
Writes stations.png / stations.pdf.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from obspy import read
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def par(key):
    for line in open("Par.file"):
        if line.split("=")[0].strip() == key:
            return float(line.split("=")[1].split("#")[0])
    raise KeyError(key)


def coords(infofile):
    """station code -> (lon, lat) from the SAC headers of the listed traces"""
    out = {}
    if not os.path.exists(infofile):
        return out
    for line in open(infofile):
        p = line.split()
        if not p:
            continue
        try:
            st = read(p[0], headonly=True)[0].stats.sac
            out[".".join(os.path.basename(p[0]).split(".")[:2])] = (st.stlo, st.stla)
        except Exception:
            pass
    return out


def ring(lon0, lat0, dist_deg, n=361):
    """points at a fixed great-circle distance from (lon0, lat0)"""
    la0, lo0, d = np.radians(lat0), np.radians(lon0), np.radians(dist_deg)
    az = np.linspace(0, 2 * np.pi, n)
    la = np.arcsin(np.sin(la0) * np.cos(d) + np.cos(la0) * np.sin(d) * np.cos(az))
    lo = lo0 + np.arctan2(np.sin(az) * np.sin(d) * np.cos(la0), np.cos(d) - np.sin(la0) * np.sin(la))
    return np.degrees(lo), np.degrees(la)


def add_base(ax, res, land=True):
    try:
        if land:      # (continent-sized fills can fail in azimuthal projections)
            ax.add_feature(cfeature.LAND.with_scale(res), facecolor="#ececec", edgecolor="none")
        ax.add_feature(cfeature.COASTLINE.with_scale(res), linewidth=0.5, edgecolor="#808080")
        if res != "110m":
            ax.add_feature(cfeature.BORDERS.with_scale(res), linewidth=0.4, edgecolor="#a0a0a0")
    except Exception:          # offline without cached Natural Earth data: plain background
        pass


evlo, evla = par("evlo"), par("evla")
P, SH = coords("stations.info"), coords("stationsSH.info")
loc = []
if os.path.exists("stationsloc.info"):
    for line in open("stationsloc.info"):
        p = line.split()
        if len(p) >= 3:
            loc.append((float(p[1]), float(p[2])))
loc = np.array(loc).reshape(-1, 2)

both = sorted(set(P) & set(SH))
ponly = sorted(set(P) - set(SH))
shonly = sorted(set(SH) - set(P))
C_BOTH, C_P, C_SH, C_LOC, C_EQ = "#3b6fb6", "#8fb3e0", "#e8853d", "#4d9a8a", "#c0392b"

fig = plt.figure(figsize=(13, 6.2))
pc = ccrs.PlateCarree()

# ---- global (teleseismic) -----------------------------------------------
axg = fig.add_axes([0.02, 0.14, 0.47, 0.76],
                   projection=ccrs.AzimuthalEquidistant(central_longitude=evlo, central_latitude=evla))
axg.set_global()
add_base(axg, "110m", land=False)
for d in (30, 60, 90):
    lo, la = ring(evlo, evla, d)
    axg.plot(lo, la, transform=ccrs.Geodetic(), color="#9a9a9a", lw=0.8, ls="--")
    lo1, la1 = ring(evlo, evla, d, n=2)
    axg.text(lo1[0], la1[0], f"{d}°", transform=pc, fontsize=9, color="#6b6b6b",
             ha="left", va="bottom")
for names, src, fc, ec, lab in ((ponly, P, "white", C_BOTH, f"P only ({len(ponly)})"),
                                 (shonly, SH, C_SH, "k", f"SH only ({len(shonly)})"),
                                 (both, P, C_BOTH, "k", f"P and SH ({len(both)})")):
    if names:
        xy = np.array([src[n] for n in names])
        axg.scatter(xy[:, 0], xy[:, 1], transform=pc, marker="^", s=50, facecolor=fc,
                    edgecolor=ec, linewidth=1.0 if fc == "white" else 0.4, label=lab, zorder=3)
axg.scatter([evlo], [evla], transform=pc, marker="*", s=260, color=C_EQ, edgecolor="k",
            linewidth=0.6, zorder=4, label="epicenter")
axg.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=9.5, frameon=False, ncol=4,
           handletextpad=0.2, columnspacing=1.0)
axg.set_title(f"Teleseismic stations: {len(P)} P, {len(SH)} SH", fontsize=12)

# ---- regional ------------------------------------------------------------
axl = fig.add_axes([0.54, 0.10, 0.43, 0.76],
                   projection=ccrs.AzimuthalEquidistant(central_longitude=evlo, central_latitude=evla))
pts = np.vstack([loc, [[evlo, evla]]]) if len(loc) else np.array([[evlo, evla]])
lo_min, lo_max = pts[:, 0].min(), pts[:, 0].max()
la_min, la_max = pts[:, 1].min(), pts[:, 1].max()
span = max(lo_max - lo_min, la_max - la_min, 2.0)
pad = 0.12 * span + 0.3
cx, cy = (lo_min + lo_max) / 2, (la_min + la_max) / 2
half = span / 2 + pad
kx = 1 / np.cos(np.radians(cy))          # equal ground distance east-west and north-south
axl.set_extent([cx - half * kx, cx + half * kx, cy - half, cy + half], crs=pc)
add_base(axl, "50m")
maxkm = max(50.0, 111.2 * np.hypot((pts[:, 0] - evlo) * np.cos(np.radians(evla)),
                                    pts[:, 1] - evla).max())
step = min([50, 100, 200, 250, 500, 1000], key=lambda v: abs(maxkm / v - 3))
for dk in np.arange(step, maxkm + step, step):
    lo, la = ring(evlo, evla, dk / 111.195)
    axl.plot(lo, la, transform=ccrs.Geodetic(), color="#9a9a9a", lw=0.7, ls="--")
    if evla + dk / 111.195 < cy + half * 0.97:      # label only rings whose top is on the map
        axl.text(evlo, evla + dk / 111.195, f"{dk:.0f} km", transform=pc, fontsize=8.5,
                 color="#6b6b6b", ha="center", va="bottom")
if len(loc):
    axl.scatter(loc[:, 0], loc[:, 1], transform=pc, marker="^", s=50, color=C_LOC,
                edgecolor="k", linewidth=0.4, zorder=3, label=f"regional ({len(loc)})")
axl.scatter([evlo], [evla], transform=pc, marker="*", s=300, color=C_EQ, edgecolor="k",
            linewidth=0.6, zorder=4, label="epicenter")
gl = axl.gridlines(draw_labels=True, linewidth=0.3, color="#c8c8c8")
gl.top_labels = gl.right_labels = False
gl.rotate_labels = False
axl.legend(loc="lower left", fontsize=9.5, framealpha=0.9)
axl.set_title(f"Regional stations: {len(loc)} (three components)", fontsize=12)

fig.savefig("stations.png", dpi=140)
fig.savefig("stations.pdf")
print(f"wrote stations.png / .pdf ({len(P)} P, {len(SH)} SH, {len(loc)} regional)")
