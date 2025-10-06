#!/usr/bin/env python3

# Credits to Giovani Gonçalves Ribamar

# This script processes metadata files from diffraction ring images.
# To run it, ensure that the script is in the same directory as:
#   - The .D5DT file
#   - The .log files
#   - A folder named "metadata" containing all metadata files
#     corresponding to the same measurement
#
# If the measurement was split into two parts, place the metadata
# files from both parts together in the same "metadata" folder.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.interpolate import interp1d
from pathlib import Path


def fetch_timestamp(fname):
    """
    Function to take the timestamp value from a Metadata file
    The metadata files comes from .tiff diffraction images
    It searches for the line containing 'timeStamp' and extracts the value
    after the '=' character
    If no such line is found, it returns None
    """
    timestamp = None
    with open(fname, 'r') as f:
        for line in f:
            if 'timeStamp' in line:
                timestamp = float(line.split('=')[-1])
                break
    return timestamp


def parse_log_files(fnames):
    """
    Parse a log file taking the informations of: (i) Name (basename);
    (ii) acquisition index (number of the image); (iii) temperature given
    by the log file (np.nan substitutes "unknown" values)
    Return a pandas DataFrame
    """
    nr1 = []
    nr2 = []
    nr = []
    mode_raw = []
    mode = []
    acq = []
    number_images = []
    for fname in fnames:
        df = pd.read_csv(fname, sep=r"\s+", comment='!', header=None)
        mode_raw += list(df[4])
        
        nr1 += list(df[7].replace('unknown', np.nan).astype(float))
        nr2 += list(df[10].replace('unknown', np.nan).astype(float))
        acq += list(df[2].replace('unknown', np.nan).astype(int))
        number_images += list(df[3].replace('unknown', np.nan).astype(int))
    for i in range(len(nr1)):
        if mode_raw[i] == "slow":
            if np.isnan(nr1[i]) and not np.isnan(nr2[i]):
                nr.append(nr2[i])
            else:
                nr.append((nr1[i] + nr2[i]) / 2)
        else:
            fast_mode_range = np.linspace(nr1[i], nr2[i], number_images[i])
            for j in range(number_images[i]):
                nr.append(fast_mode_range[j])
    return pd.DataFrame(dict(nr=nr))


if __name__ == '__main__':
    print("Creating a synchronized file between diffraction metadata and dilatometry files...")
    print(".\n.\n.")
    log_files = Path().glob("*.log")

    # Fetches metadata from all points, including fast acquisition
    timestamp_metadata = OrderedDict()
    for fpath in sorted(Path('metadata').glob('*.tif.metadata')):
        timestamp = fetch_timestamp(fpath)
        if timestamp is None:
            print(f'Could not find timeStamp in file {str(fpath)}')
            continue
        timestamp_metadata[str(fpath)] = timestamp

    
    # Reads dilatometry data
    dil_file = list(Path().glob("*.D5DT"))[0]

    df_dil = pd.read_csv(dil_file,
                         header=None, sep=r"\s+")
    time2temp = interp1d(df_dil[0], df_dil[1], bounds_error=False,
                         fill_value=np.nan)

    print("Dilatometer files read!")
    
    
    print("Synchronizing...")
    print(".\n.\n.")
    df_temp = parse_log_files(log_files)
    df_time_values = []
    df_temp_values = []
    df_change_in_lenght = []
    df_force = []
    
    acqindex2time = interp1d(df_dil.index, df_dil[0],
                    bounds_error=False, fill_value=np.nan)
    acqindex2temp = interp1d(df_dil.index, df_dil[1],
                    bounds_error=False, fill_value=np.nan)
    acqindex2change_in_lenght = interp1d(df_dil.index, df_dil[2],
                    bounds_error=False, fill_value=np.nan)
    acqindex2force = interp1d(df_dil.index, df_dil[3],
                    bounds_error=False, fill_value=np.nan)
    
    df_temp["time"] = acqindex2time(df_temp["nr"])
    df_temp["temp"] = acqindex2temp(df_temp["nr"])
    df_temp["change_in_lenght"] = acqindex2change_in_lenght(df_temp["nr"])
    df_temp["force"] = acqindex2force(df_temp["nr"])

    
    print("Wait for the plots")
    print(".\n.\n.")
    # Plots stuff
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    
    ax1.plot(df_dil[0], df_dil[1], 'r-', label="Dilatometry data")
    ax1.plot(df_temp['time'], df_temp['temp'], 'kx', label="Diffraction data")
    ax1.set_xlabel("Time (s)", size=15)
    ax1.set_ylabel("Temperature (°C)", size=15)
    ax1.legend()
    
    ax2.plot(df_dil[0], df_dil[2], 'r-', label="Dilatometry data")
    ax2.plot(df_temp['time'], df_temp['change_in_lenght'], 'kx', label="Diffraction data")
    ax2.set_xlabel("Time (s)", size=15)
    ax2.set_ylabel(r"Change in Lenght ($\mu$m)", size=15)
    ax2.legend()
    
    ax3.plot(df_dil[0], df_dil[3], 'r-', label="Dilatometry data")
    ax3.plot(df_temp['time'], df_temp['force'], 'kx', label="Diffraction data")
    ax3.set_xlabel("Time (s)", size=15)
    ax3.set_ylabel("Force (N)", size=15)
    ax3.legend()

    sync_file = pd.DataFrame()
    sync_file["Image"] = [k.lstrip("metadata\\").rstrip(".tif.metadata")
                          for k in timestamp_metadata.keys()]
    sync_file["Time(s)"] = df_temp["time"]
    sync_file["Temperature(oC)"] = df_temp["temp"]
    sync_file["Change_in_Length(micrometer)"] = df_temp["change_in_lenght"]
    sync_file["Force(N)"] = df_temp["force"]    

    name = sync_file["Image"][0].rstrip("_1-00001")

    sync_file.to_csv("./" + name + "_sync_file.txt", sep=";", index=None, encoding="utf-8")

    fig.tight_layout()
    fig.savefig("./" + name + "_plot_sync_results.jpg", dpi=400)
    
    print("Just close the figure window to end the program.")
    print(".\n.\n.")
    print("Your created synchronized file is in the same folder of this executable file")
    print("Bye!!!")
    
    plt.show()
