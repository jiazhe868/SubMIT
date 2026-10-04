import os
import sys
import requests
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader

# Function to download NEIC earthquake aftershock catalog
def download_aftershock_catalog(mainshock_time, longitude, latitude, depth, mw=None):
    try:
        mainshock_time_parsed = datetime.strptime(mainshock_time, "%Y-%m-%dT%H%M%S")
    except ValueError:
    # If the first format fails, try the second format
        try:
            mainshock_time_parsed = datetime.strptime(mainshock_time, "%Y/%m/%dT%H:%M:%S.%f")
        except ValueError:
        # Handle the case where both formats fail
            print("Time format not recognized.")
            mainshock_time_parsed = None
    mainshock_time = mainshock_time_parsed
    end_time = mainshock_time + timedelta(days=7)
    
    # Define search box: scale with the Wells-Coppersmith rupture length so
    # giant events keep their full aftershock zone (M8.8 -> ~850 km half-width;
    # the old fixed 200 km clipped Kamchatka's zone at 260 km); floor 200 km.
    # caldens.py still handles outliers and sparse catalogs downstream.
    half_km = max(200.0, 1.5 * 10 ** (0.59 * mw - 2.44)) if mw else 200.0
    min_latitude = latitude - (half_km / 111)
    max_latitude = latitude + (half_km / 111)
    min_longitude = longitude - (half_km / (111 * abs(np.cos(np.radians(latitude)))))
    max_longitude = longitude + (half_km / (111 * abs(np.cos(np.radians(latitude)))))
    min_depth = max(depth - 30, 0)  # 30 km shallower, ensure depth is not negative
    max_depth = depth + 30  # 30 km deeper

    # Prepare the query URL
    url = (
        "https://earthquake.usgs.gov/fdsnws/event/1/query"
        "?format=geojson"
        f"&starttime={mainshock_time.isoformat()}"
        f"&endtime={end_time.isoformat()}"
        f"&minlatitude={min_latitude}"
        f"&maxlatitude={max_latitude}"
        f"&minlongitude={min_longitude}"
        f"&maxlongitude={max_longitude}"
        f"&mindepth={min_depth}"
        f"&maxdepth={max_depth}"
        f"&minmagnitude=0"  # Adjust minimum magnitude if needed
    )

    # Make the request to USGS API; on failure, leave an empty catalog rather
    # than aborting the whole preparation pipeline (caldens.py falls back to a
    # magnitude-scaled uniform prior box)
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        data = response.json()
    except Exception as exc:
        print(f"WARNING: aftershock query failed ({exc}); writing empty catalog.")
        data = {'features': []}

    # Extract earthquake details
    earthquakes = []
    for feature in data['features']:
        properties = feature['properties']
        geometry = feature['geometry']
        time = datetime.utcfromtimestamp(properties['time'] / 1000).strftime("%Y-%m-%dT%H:%M:%S")
        mag = properties['mag']
        lon, lat, dep = geometry['coordinates']

        # Skip the mainshock itself and duplicate agency listings of it:
        # within 120 s of origin, within 25 km, and not clearly smaller
        t_ev = datetime.strptime(time, "%Y-%m-%dT%H:%M:%S")
        d_km = np.hypot((lon - longitude) * 111.32 * np.cos(np.radians(latitude)),
                        (lat - latitude) * 110.574)
        if abs((t_ev - mainshock_time).total_seconds()) < 120 and d_km < 25 \
                and (mw is None or (mag is not None and mag >= mw - 0.5)):
            continue

        earthquakes.append([lon, lat, dep, mag, time])

    # Convert to DataFrame
    df = pd.DataFrame(earthquakes, columns=['Longitude', 'Latitude', 'Depth', 'Magnitude', 'Origin Time'])

    # Save to an ASCII file
    df.to_csv('catalog.txt', sep=' ', index=False, header=False, float_format='%.4f')
    print("Aftershock catalog saved to 'catalog.txt'")
    # Save only the Longitude and Latitude columns to another ASCII file.
    # Magnitude floor min(4.0, max(2.5, Mw-4)): footprint-defining aftershocks
    # of an M7 are M>~3; permanent background microseismicity otherwise
    # dominates the catalog (2024 Mendocino M7.0: 47% of the raw box catalog
    # was Geysers geothermal M<2 swarm 250 km SE -> inflated the median so
    # caldens outlier rejection could not fire, and falsely triggered the
    # extended-source branch). Cap 4.0 keeps M8+ catalogs (M>=4.1) intact.
    dffp = df
    if mw is not None:
        mfloor = min(4.0, max(2.5, float(mw) - 4.0))
        dffp = df[df['Magnitude'] >= mfloor]
        ndrop = len(df) - len(dffp)
        if ndrop:
            print(f"fetch_aftershocks: magnitude floor M>={mfloor:.1f} drops "
                  f"{ndrop} of {len(df)} events (background microseismicity)")
    dffp[['Longitude', 'Latitude']].to_csv('lola.dat', sep=' ', index=False, header=False, float_format='%.4f')

    # BACKGROUND catalog (same box/depth/magnitude floor, 365 d ending 1 d
    # before the mainshock): caldens.py keeps an aftershock only where the
    # 7-day rate significantly exceeds this background rate (Poisson test) -
    # permanent clusters (Geysers geothermal swarm, triple-junction background)
    # cancel out regardless of magnitude. On fetch failure no lola_bg.dat is
    # written and caldens falls back to the unfiltered behavior.
    try:
        bg_end = mainshock_time - timedelta(days=1)
        bg_start = bg_end - timedelta(days=365)
        bg_url = (
            "https://earthquake.usgs.gov/fdsnws/event/1/query"
            "?format=geojson"
            f"&starttime={bg_start.isoformat()}"
            f"&endtime={bg_end.isoformat()}"
            f"&minlatitude={min_latitude}"
            f"&maxlatitude={max_latitude}"
            f"&minlongitude={min_longitude}"
            f"&maxlongitude={max_longitude}"
            f"&mindepth={min_depth}"
            f"&maxdepth={max_depth}"
            "&minmagnitude=0"
        )
        bg_resp = requests.get(bg_url, timeout=60)
        bg_resp.raise_for_status()
        bg_events = [[f['geometry']['coordinates'][0], f['geometry']['coordinates'][1],
                      f['geometry']['coordinates'][2], f['properties']['mag']]
                     for f in bg_resp.json()['features'] if f['properties']['mag'] is not None]
        bg_df = pd.DataFrame(bg_events, columns=['Longitude', 'Latitude', 'Depth', 'Magnitude'])
        if mw is not None:
            bg_df = bg_df[bg_df['Magnitude'] >= min(4.0, max(2.5, float(mw) - 4.0))]
        bg_df[['Longitude', 'Latitude']].to_csv('lola_bg.dat', sep=' ', index=False,
                                                header=False, float_format='%.4f')
        print(f"fetch_aftershocks: background catalog {len(bg_df)} events/365 d "
              f"(same box + magnitude floor) -> lola_bg.dat")
    except Exception as exc:
        print(f"WARNING: background-catalog fetch failed ({exc}); "
              "rate filtering disabled for this run.")

    # Mainshock reference for downstream scripts (caldens.py reads the Mw)
    with open('mainshock.dat', 'w') as f:
        f.write(f"{longitude:.4f} {latitude:.4f} {depth:.2f} {mw if mw is not None else 'nan'}\n")

    # Plot aftershocks and mainshock (optional; never fatal for the pipeline)
    try:
        plot_aftershocks(df, longitude, latitude, min_longitude, max_longitude, min_latitude, max_latitude)
    except Exception as exc:
        print(f"WARNING: aftershock map plotting failed ({exc}); continuing.")

# Function to plot the aftershocks and mainshock
def plot_aftershocks(df, mainshock_lon, mainshock_lat, min_lon, max_lon, min_lat, max_lat):
    # Create a map with topography, bathymetry, and coastlines
    fig = plt.figure(figsize=(12, 10))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Add topography and bathymetry
    ax.stock_img()  # Provides a basic topography and bathymetry background
    
    # Add coastlines
    ax.add_feature(cfeature.COASTLINE, linewidth=1)
    
    # Plot the aftershocks
    sc = ax.scatter(
        df['Longitude'], df['Latitude'], 
        s=df['Magnitude']**2 * 5,  # Marker size proportional to magnitude squared
        c=df['Depth'], cmap='viridis', 
        edgecolor='k', label='Aftershocks', transform=ccrs.PlateCarree()
    )

    # Plot the mainshock as a star
    ax.plot(mainshock_lon, mainshock_lat, 'r*', markersize=15, label='Mainshock', transform=ccrs.PlateCarree())

    # Set plot limits to match the search area
    ax.set_extent([min_lon, max_lon, min_lat, max_lat])
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = plt.MultipleLocator(0.5)
    gl.ylocator = plt.MultipleLocator(0.5)
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

    # Add colorbar and labels
    cbar = plt.colorbar(sc, ax=ax, orientation='vertical', shrink=0.5, pad=0.02)
    cbar.set_label('Depth (km)')
    ax.set_title('Aftershock Map with Mainshock and Topography')
    ax.legend()

    # Save the figure
    plt.savefig('aftershock_map.pdf', format='pdf')
    #plt.show()

# Main execution
if __name__ == "__main__":
    # Extract arguments
    mainshock_time = sys.argv[1]  # Example: "2024-08-19T123000"
    longitude = float(sys.argv[2])  # Example: 121.5837
    latitude = float(sys.argv[3])  # Example: 23.8607
    depth = float(sys.argv[4])  # Example: 10.0
    mw = float(sys.argv[5]) if len(sys.argv) > 5 else None  # mainshock Mw (optional)

    # SUBMIT_CATALOG_CACHE=<dir>: replay a saved catalog instead of querying
    # USGS (exact reproduction of a published run; the live catalog is revised
    # over time, which shifts the aftershock prior)
    cache = os.environ.get("SUBMIT_CATALOG_CACHE")
    if cache and os.path.exists(os.path.join(cache, "lola.dat")):
        import shutil
        for f in ("catalog.txt", "lola.dat", "lola_bg.dat", "mainshock.dat"):
            if os.path.exists(os.path.join(cache, f)):
                shutil.copy(os.path.join(cache, f), f)
        print(f"fetch_aftershocks: replayed cached catalog from {cache}")
        sys.exit(0)

    # Download the aftershock catalog
    download_aftershock_catalog(mainshock_time, longitude, latitude, depth, mw)

