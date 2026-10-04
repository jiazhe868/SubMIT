import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Generate a smoothed density map of earthquakes.")
parser.add_argument('evlo', type=float, help='Hypocenter longitude')
parser.add_argument('evla', type=float, help='Hypocenter latitude')
parser.add_argument('--floor', type=float, default=0.0,
                    help='uniform density floor (0-1) added inside the grid; '
                         'useful for sparse catalogs (default 0 = original behavior)')
parser.add_argument('--nmin', type=int, default=10,
                    help='minimum usable aftershocks for a density prior; below this '
                         'a uniform bounded prior is used instead (default 10)')
parser.add_argument('--half', type=float, default=100.0,
                    help='half-width (km) of the uniform prior box used when '
                         'NO usable aftershocks exist (default 100)')
parser.add_argument('--buffer', type=float, default=30.0,
                    help='buffer (km) added around the aftershock extent when few '
                         'aftershocks define the uniform search box (default 30)')

args = parser.parse_args()

# Assign evlo and evla from arguments
evlo = args.evlo
evla = args.evla

# Load earthquake data (robust to 0- or 1-event catalogs)
try:
    seislola = np.loadtxt('lola.dat', ndmin=2)
except Exception:
    seislola = np.zeros((0, 2))
if seislola.size == 0:
    seislola = np.zeros((0, 2))

delta = 4  # grid interval in km

# Convert coordinates (hypocenter at the origin)
seisxy = seislola.copy()
if len(seislola) > 0:
    seisxy[:, 0] = (seislola[:, 0] - evlo) * 111.32 * np.cos(np.radians(evla))
    seisxy[:, 1] = (seislola[:, 1] - evla) * 110.574

# Background-rate filter (needs lola_bg.dat from fetch_aftershocks.py): keep a
# candidate aftershock only where the 7-day post-event count within R=20 km is
# a statistically significant (Poisson p<0.01) excess over the local 365-day
# background rate. Permanent sources - geothermal swarms (2024 Mendocino: the
# Geysers field 250 km SE), triple-junction background - run at their steady
# rate during the aftershock week and cancel; true aftershock zones run at
# 10-100x background and pass. Isolated events with zero local background pass
# (the median outlier rejection below still governs them).
if len(seisxy) > 0:
    try:
        bglola = np.loadtxt('lola_bg.dat', ndmin=2)
    except Exception:
        bglola = np.zeros((0, 2))
    if bglola.size > 0:
        from scipy.stats import poisson
        bgxy = bglola.copy()
        bgxy[:, 0] = (bglola[:, 0] - evlo) * 111.32 * np.cos(np.radians(evla))
        bgxy[:, 1] = (bglola[:, 1] - evla) * 110.574
        R = 20.0
        keep = np.ones(len(seisxy), dtype=bool)
        for i in range(len(seisxy)):
            a = int(np.sum(np.hypot(seisxy[:, 0] - seisxy[i, 0],
                                    seisxy[:, 1] - seisxy[i, 1]) <= R))
            nbg = int(np.sum(np.hypot(bgxy[:, 0] - seisxy[i, 0],
                                      bgxy[:, 1] - seisxy[i, 1]) <= R))
            expected = nbg * 7.0 / 365.0
            # P(X >= a | expected) < 0.01 -> significant activation
            keep[i] = poisson.sf(a - 1, expected) < 0.01 if expected > 0 else True
        ndrop = int((~keep).sum())
        if ndrop:
            print(f"caldens: background-rate filter dropped {ndrop} of "
                  f"{len(seisxy)} candidate aftershocks (steady-rate sources)")
        seisxy = seisxy[keep]

# Reject outliers data-driven (no rupture-length scaling: true rupture size is
# unknown a priori - that is what the inversion determines). An aftershock is an
# outlier if its distance from the hypocenter exceeds
# max(median + 3*1.4826*MAD, 50 km) of the distance distribution.
if len(seisxy) > 0:
    dist_hyp = np.hypot(seisxy[:, 0], seisxy[:, 1])
    med = np.median(dist_hyp)
    mad = np.median(np.abs(dist_hyp - med))
    rmax = max(med + 3 * 1.4826 * mad, 50.0)
    nrej = int(np.sum(dist_hyp > rmax))
    if nrej:
        print(f"caldens: rejected {nrej} outlier aftershock(s) beyond {rmax:.0f} km "
              f"(median dist {med:.0f} km)")
    seisxy = seisxy[dist_hyp <= rmax]

if len(seisxy) <= args.nmin:
    # <= nmin aftershocks: a density estimated from a handful of events is not
    # informative -> UNIFORM prior. But the aftershock LOCATIONS still carry
    # search-boundary information: box = aftershock extent + --buffer km,
    # always also covering +/-21 km around the hypocenter. With zero usable
    # aftershocks, fall back to +/---half km around the hypocenter.
    if len(seisxy) > 0:
        buf = args.buffer
        minx = np.floor((min(np.min(seisxy[:, 0]) - buf, -21.0)) / delta) * delta - 1
        maxx = np.ceil((max(np.max(seisxy[:, 0]) + buf, 21.0)) / delta) * delta + 1
        miny = np.floor((min(np.min(seisxy[:, 1]) - buf, -21.0)) / delta) * delta - 1
        maxy = np.ceil((max(np.max(seisxy[:, 1]) + buf, 21.0)) / delta) * delta + 1
        print(f"caldens: only {len(seisxy)} usable aftershock(s) (<= nmin={args.nmin}); "
              f"uniform prior over aftershock extent +/-{buf:.0f} km buffer: "
              f"x [{minx:.0f},{maxx:.0f}], y [{miny:.0f},{maxy:.0f}]")
    else:
        half = args.half
        print(f"caldens: no usable aftershocks; uniform prior within "
              f"+/-{half:.0f} km of the hypocenter")
        minx = -np.ceil(half / delta) * delta - 1
        maxx = np.ceil(half / delta) * delta + 1
        miny, maxy = minx, maxx
else:
    minx = np.round(np.min(seisxy[:, 0]) / delta) * delta - 21  # -/+20 km extension
    maxx = np.round(np.max(seisxy[:, 0]) / delta) * delta + 21
    miny = np.round(np.min(seisxy[:, 1]) / delta) * delta - 21
    maxy = np.round(np.max(seisxy[:, 1]) / delta) * delta + 21
    # the search region must contain the hypocenter (x = y = 0) comfortably
    minx = min(minx, -21.0); maxx = max(maxx, 21.0)
    miny = min(miny, -21.0); maxy = max(maxy, 21.0)

nx = round((maxx - minx) / delta)
ny = round((maxy - miny) / delta)

xx = seisxy[:, 0]
yy = seisxy[:, 1]

# Filter earthquakes within the desired area
index = (xx > minx) & (xx < maxx) & (yy > miny) & (yy < maxy)
seisxy = seisxy[index, :]
xx = xx[index]
yy = yy[index]

# Compute histogram
n, xedges, yedges = np.histogram2d(xx, yy, bins=[np.arange(minx, maxx + delta, delta), np.arange(miny, maxy + delta, delta)])
n1 = n.T

# Smooth the density map
n2 = gaussian_filter(n1, sigma=2)

# Normalize the smoothed density
maxn = np.max(n2)
# If <= nmin usable aftershocks, use a flat prior over the fallback box
if len(seisxy) <= args.nmin or maxn <= 0:
    n2[:] = 1  # Set all values in n2 to 1
else:
    # Apply thresholds to n2 based on maxn
    n2[n2 > maxn / 4] = maxn / 4
    # CONNECTED corridor: sparse, separated aftershock clusters (e.g. a doublet
    # with a quiet gap between the two rupture zones) would otherwise be zeroed
    # in between, and the x/y random walk could never traverse east-west.
    # Bridge the catalog along minimum-spanning-tree edges and keep a low
    # density floor along that corridor; hard-zero only OFF-corridor bins.
    bridge = np.zeros_like(n2)
    try:
        from scipy.sparse.csgraph import minimum_spanning_tree
        from scipy.spatial.distance import cdist
        from scipy.ndimage import label as nd_label, distance_transform_edt
        pts = np.unique(np.round(np.c_[xx, yy] / delta) * delta, axis=0)
        if len(pts) > 1:
            # bridges exist to CONNECT separated aftershock clusters (the
            # Venezuela doublet case: depleted gap between two rupture zones),
            # NOT to skeletonize the interior of one cluster or to grant
            # aftershock-zone probability to no-aftershock locations. Keep
            # only MST edges whose endpoints lie in DIFFERENT connected
            # components of the density support - within a cluster the
            # heatmap itself already provides the path.
            # label on the SURVIVING support (post final threshold
            # maxn/20): smoothing tails can weakly connect two real
            # clusters at n2>0 level, which would suppress the needed
            # bridge and then be cut by the threshold -> disconnected
            comp, ncomp = nd_label(n2 > maxn / 20.0)
            # keep only MAJOR components (>=5% of total density mass): a giant
            # rupture's aftershock box sweeps in permanent subduction-zone
            # seismicity that survives the magnitude/rate filters as scattered
            # speckle; bridging every speckle rebuilt the spiky tree
            # (Kamchatka: 93 bridges) and let subevents sit on unsupported
            # spots. Minor components are ZEROED, not bridged.
            if ncomp > 1:
                mass = np.array([n2[comp == k].sum() for k in range(1, ncomp + 1)])
                minor = [k + 1 for k in range(ncomp) if mass[k] < 0.05 * mass.sum()]
                if minor:
                    for k in minor:
                        n2[comp == k] = 0.0
                    comp, ncomp = nd_label(n2 > maxn / 20.0)
                    print(f"caldens: zeroed {len(minor)} minor density component(s) "
                          f"(<5% mass; background speckle), {ncomp} major kept")
            plab = np.zeros(len(pts), dtype=int)
            for k in range(len(pts)):
                ix = min(max(int((pts[k, 0] - minx) / delta), 0), n2.shape[1] - 1)
                iy = min(max(int((pts[k, 1] - miny) / delta), 0), n2.shape[0] - 1)
                plab[k] = comp[iy, ix]
            T = minimum_spanning_tree(cdist(pts, pts)).tocoo()
            nbridge = 0
            for i, j in zip(T.row, T.col):
                if plab[i] == plab[j] and plab[i] != 0:
                    continue
                pq = pts[j] - pts[i]
                nseg = max(2, int(np.hypot(*pq) / (delta / 2.0)))
                for t in np.linspace(0.0, 1.0, nseg):
                    x0, y0 = pts[i] + t * pq
                    ix = int((x0 - minx) / delta)
                    iy = int((y0 - miny) / delta)
                    if 0 <= iy < n2.shape[0] and 0 <= ix < n2.shape[1]:
                        bridge[iy, ix] = 1.0
                nbridge += 1
            if nbridge:
                # smooth half-strength ZONE, not a line: Gaussian falloff from
                # the connecting path, peak 0.5 of the clipped cluster level,
                # sigma 3 cells (~12 km at delta=4), cut at 0.1. Traversable by
                # the random walk (prior ratio cluster->bridge = 0.5) but a
                # subevent PARKED on a bridge pays a persistent prior deficit.
                dcells = distance_transform_edt(bridge == 0)
                bridge = 0.5 * np.exp(-0.5 * (dcells / 3.0) ** 2)
                bridge[bridge < 0.1] = 0.0
                print(f"caldens: {nbridge} inter-cluster bridge(s), half-strength smooth zones")
            else:
                bridge[:] = 0.0
    except Exception as exc:
        print(f"caldens: WARNING corridor bridge failed ({exc})")
    # combine: heatmap everywhere, half-strength smooth bridges only between
    # separated clusters; off-support and off-bridge stays zero
    n2 = np.maximum(n2, (maxn / 4) * bridge)
    n2[(n2 < maxn / 20) & (bridge <= 0.0)] = 0
    # optional uniform floor so sparse catalogs do not zero out the search area
    if args.floor > 0:
        n2 = np.maximum(n2, args.floor * np.max(n2))

# Prepare grid for final output
A, B = np.meshgrid(xedges[:-1], yedges[:-1])
final = np.vstack([A.ravel(), B.ravel(), n2.ravel() / np.max(n2)]).T
np.savetxt('seisdens.dat', final, fmt='%.2f', delimiter=' ')  # Smoothed aftershock density
np.savetxt('seisgrids.dat', [minx, maxx, miny, maxy, delta], fmt='%.2f', delimiter=' ')  # Grid parameters

# Determine the min and max x, y where n2 > 0 (guard: never an empty region)
nonzero_indices = np.where(n2 > 0)
if len(nonzero_indices[0]) == 0:
    print("caldens: WARNING - density all zero, using full grid as edges")
    n2[:] = 1
    nonzero_indices = np.where(n2 > 0)
minx, maxx = xedges[nonzero_indices[1]].min(), xedges[nonzero_indices[1]].max()
miny, maxy = yedges[nonzero_indices[0]].min(), yedges[nonzero_indices[0]].max()

# Save the min and max x, y values to 'edges.dat'
with open('edges.dat', 'w') as f:
    f.write(f"{minx:.2f} {maxx:.2f} {miny:.2f} {maxy:.2f}")

# Plot the density map and save as PDF
plt.figure()
plt.pcolor(xedges[:-1], yedges[:-1], n2, shading='auto')
plt.colorbar()
plt.plot(0, 0, 'p', markersize=10, markerfacecolor='red')
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Aftershock Density in Map View')
plt.savefig('aftershock_density.pdf', format='pdf')
#plt.show()
