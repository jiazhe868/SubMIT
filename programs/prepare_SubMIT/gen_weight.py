import sys
import numpy as np
from obspy import read

# Function to calculate duration based on magnitude Mw
def calculate_duration(mw):
    return 5 * 10 ** (0.5 * (mw - 6))

# Function to calculate the window for P and SH waves
def calculate_window_p_sh(t1_index, sampling_rate, duration):
    start_index = int(t1_index - 20 * sampling_rate)  # 20 seconds before t1
    end_index = int(t1_index + (duration + 70) * sampling_rate)  # Duration + 70 seconds after t1
    return start_index, end_index

# Function to calculate the window for regional waves
def calculate_window_regional(t1_index, sampling_rate):
    start_index = int(t1_index - 20 * sampling_rate)  # 20 seconds before t1
    end_index = int(t1_index + 200 * sampling_rate)  # 200 seconds after t1
    return start_index, end_index

# Function to read and process SAC files, apply bandpass filter, and calculate summarized L2 norm
def calculate_summarized_l2_norm(file_list, window_func, sampling_rate, duration=None, freqmin=0.02, freqmax=0.2):
    combined_data = []
    for sac_file in file_list:
        try:
            st = read(sac_file)  # Read SAC file
            tr = st[0]  # Get the trace
            tr.detrend(type='demean')  # Remove mean
            tr.detrend(type='linear')  # Remove linear trend
            tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=True)  # Apply bandpass filter

            t1 = tr.stats.sac.t1  # Get the t1 marker time in seconds
            npts = tr.stats.npts  # Number of points in the data
            t1_index = int((t1 - tr.stats.sac.b) * sampling_rate)  # Convert t1 to sample index

            # Determine window based on the function passed
            if duration:
                start_index, end_index = window_func(t1_index, sampling_rate, duration)
            else:
                start_index, end_index = window_func(t1_index, sampling_rate)
            
            # Ensure the indices are within the bounds of the data
            start_index = max(0, start_index)
            end_index = min(npts, end_index)

            # Trim the data using calculated indices
            trimmed_data = tr.data[start_index:end_index]
            
            # Append the trimmed data to the combined list
            combined_data.extend(trimmed_data)
        except Exception as e:
            print(f"Error processing {sac_file}: {e}")
    
    # Calculate summarized L2 norm for the combined data
    combined_data = np.array(combined_data)  # Convert combined data to a numpy array
    summarized_l2_norm = np.linalg.norm(combined_data)  # Calculate summarized L2 norm
    return summarized_l2_norm

# Function to read filenames from info files
def read_filenames(file_path):
    filenames = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines
                filenames.append(line.split()[0])  # Read first column (SAC filenames)
    return filenames

# Function to calculate weights with two significant digits
def calculate_weights(l2_norms):
    weights = 1 / np.array(l2_norms)  # Inverse of L2 norms
    weights /= weights[1]  # Normalize SH weight to 1 (SH is assumed to be at index 1)
    # Regional trust factor 1/3: the 1D velocity model is least reliable for
    # regional full waveforms. Calibrated against BOTH production hand-tunings:
    # Venezuela band-corrected 0.32 -> 0.11 (manual 0.1), Calama/Myanmar-style
    # 0.056 -> 0.019 (manual 0.02).
    weights[2] /= 3.0

    
    # Function to format to two significant digits
    def format_two_sig_digits(x):
        if x == 0:
            return 0.0
        else:
            return float(f"{x:.2g}")

    # Apply formatting to weights
    weights = [format_two_sig_digits(w) for w in weights]
    
    return weights

# Main script
if __name__ == "__main__":
    # Ensure magnitude Mw is passed as an argument
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <Mw>")
        sys.exit(1)

    try:
        mw = float(sys.argv[1])  # Convert argument to float
    except ValueError:
        print("Magnitude Mw should be a floating point number.")
        sys.exit(1)
    
    # Calculate duration based on magnitude Mw
    duration = calculate_duration(mw)
    
    # Paths to info files
    p_wave_info = 'stations.info'
    sh_wave_info = 'stationsSH.info'
    regional_wave_info = 'stationsloc.info'

    # Read SAC filenames from info files
    p_wave_files = read_filenames(p_wave_info)
    sh_wave_files = read_filenames(sh_wave_info)

    # For regional waves, append 'e', 'n', 'z' components to base filenames
    regional_wave_files = []
    base_names = read_filenames(regional_wave_info)
    for base_name in base_names:
        for comp in ['e', 'n', 'z']:
            regional_wave_files.append(f"{base_name}{comp}")

    # Example to get sampling rate from the first file
    example_trace = read(p_wave_files[0])[0]
    sampling_rate = example_trace.stats.sampling_rate

    # Per-type norms in each type's ACTUAL inversion band (same formulas as
    # gen_par.py - keep in sync). The old hardcoded 0.02-0.2 Hz ignored the
    # filtering the inversion applies: P/SH include 0.005-0.02 Hz where SH
    # carries relatively more energy, so SH was effectively overweighted and
    # could drag mechanisms away from the P radiation pattern.
    hf_body = 0.2
    # aperture-limited Rayleigh corner - same piecewise rule as gen_par.py
    # (hypocenter from mainshock.dat: lon lat dep mw)
    import math as _m
    max_sta = 0.0
    try:
        _p = open('mainshock.dat').read().split()
        evlo, evla = float(_p[0]), float(_p[1])
        for line in open('stationsloc.info'):
            t = line.split()
            if len(t) >= 3:
                _dx = (float(t[1]) - evlo) * 111.32 * _m.cos(_m.radians(evla))
                _dy = (float(t[2]) - evla) * 110.574
                max_sta = max(max_sta, _m.hypot(_dx, _dy))
    except Exception:
        pass
    if max_sta <= 0 or max_sta <= 450:
        hf_rayl = 0.15
    elif max_sta <= 570:
        hf_rayl = round(0.15 - 0.02 * (max_sta - 450) / 120.0, 3)
    else:
        # steep: 0.069 already halves the MT amplitudes on 731-km paths
        hf_rayl = round(max(0.05, 0.13 - 0.08 * (max_sta - 570) / 130.0), 3)
    p_wave_summarized_l2_norm = calculate_summarized_l2_norm(p_wave_files, calculate_window_p_sh, sampling_rate, duration, freqmin=0.005, freqmax=hf_body)
    sh_wave_summarized_l2_norm = calculate_summarized_l2_norm(sh_wave_files, calculate_window_p_sh, sampling_rate, duration, freqmin=0.005, freqmax=hf_body)
    regional_wave_summarized_l2_norm = calculate_summarized_l2_norm(regional_wave_files, calculate_window_regional, sampling_rate, freqmin=0.02, freqmax=hf_rayl)

    # Calculate weights inversely proportional to L2 norms
    l2_norms = [p_wave_summarized_l2_norm, sh_wave_summarized_l2_norm, regional_wave_summarized_l2_norm]
    weights = calculate_weights(l2_norms)
    # a regional set of <= 2 stations cannot constrain anything (Kamchatka
    # M8.8: one filler station kept only because sub_forward requires
    # weightRayl > 0) - keep it minimally positive
    if len(base_names) <= 2:
        weights[2] = min(weights[2], 0.002)
        print(f"gen_weight: only {len(base_names)} regional station(s) - "
              f"weightRayl set to {weights[2]} (filler mode)")

    # Output weights
    print("Weights (P:SH:Regional):", weights[0], ":", weights[1], ":", weights[2])

    # Save weights to ASCII file
    with open('weights.dat', 'w') as f:
        f.write(f"{weights[0]} {weights[1]} {weights[2]}\n")

