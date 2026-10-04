import scipy.signal
if not hasattr(scipy.signal, 'hann'):
    import scipy.signal.windows
    scipy.signal.hann = scipy.signal.windows.hann

import numpy as np
import obspy
from obspy.core import UTCDateTime
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy.io.sac import SACTrace

# Load the event location, depth, and origin time
with open('loladep.dat', 'r') as file:
    event_info = file.readline().strip().split()
    event_longitude = float(event_info[0])
    event_latitude = float(event_info[1])
    event_depth = float(event_info[2])  # in kilometers
    event_origin_time = UTCDateTime(event_info[3].replace('T', '').replace('Z', ''))  # parse the time

# Read the list of data files
with open('data_list.dat', 'r') as file:
    data_files = [line.strip() for line in file.readlines()]

# Function to process each file
def process_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract header information
    start_time = UTCDateTime(lines[0].split()[4])  # beginning of the time trace
    samples_per_second = float(lines[1].split()[4])
    total_data_points = int(lines[2].split()[5])
    station_name = lines[3].split()[2]
    channel_name = lines[3].split()[4]
    station_latitude = float(lines[4].split()[2])
    station_longitude = float(lines[4].split()[4])
    dt = 1/samples_per_second
    time_difference = start_time - event_origin_time

    # Extract data
    data = np.array([float(line.strip()) for line in lines[6:]])
    
    # Create a Trace
    trace = Trace(data=data)
    trace.stats.starttime = start_time
    trace.stats.station = station_name
    trace.stats.channel = channel_name
    trace.stats.delta = dt
    trace.integrate(method='cumtrapz')
    trace.detrend(type='linear')  # Remove linear trend
    trace.detrend(type='demean')  # Remove mean
    trace.filter('highpass', freq=0.01, corners=4, zerophase=True)  # High-pass filter with zero-phase
    trace.taper(max_percentage=0.05, type='hann')  # Taper the first and last 5% of the data
    trace.decimate(factor=10)  # Decimate the data by a factor of 100

    sac_trace = SACTrace.from_obspy_trace(trace)
    sac_trace.evlo = event_longitude
    sac_trace.evla = event_latitude
    sac_trace.stlo = station_longitude
    sac_trace.stla = station_latitude
    sac_trace.evdp = event_depth
    sac_trace.stdp = 0
    sac_trace.lcalda = True
    sac_trace.o = 0
    sac_trace.b = time_difference
    sac_trace.t1 = 0
   
    component_code = 'z' if channel_name.endswith('Z') else 'n' if channel_name.endswith('N') else 'e'
    sacname = file_path.replace('.txt', '.sac')
    sacname = f"C1.{station_name}..{component_code}"

    sac_trace.write(sacname)
    print('finishing: ',sacname)

# Process each file
for data_file in data_files:
    process_file(data_file)

