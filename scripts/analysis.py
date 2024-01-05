import argparse
import glob

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

import os

def create_directory(directory_path):
    # Check if the directory exists
    if not os.path.exists(directory_path):
        # If it doesn't exist, create the directory
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")



def examine_stats(quad_clean):
    pass



orbit = 43416
binning = 10
destination_path = "../../data/evt_files"
pattern_evt = f"*{orbit}*_quad_clean.evt"
plot_path = "../plots"
evt_files = glob.glob(f'{destination_path}/*{pattern_evt}')

create_directory(f"{plot_path}/{orbit}")

time_stamps = []
evt_file = evt_files[0]
with fits.open(evt_file) as hdul:
    for quarter in range(1,5):
        q_data = hdul[quarter].data['Time']
        start = hdul[quarter].header['TSTARTI']
        end = hdul[quarter].header['TSTOPI']
        q_data = q_data[(q_data>start) & (q_data<end)]
        time_stamps.append(q_data)
        print(start,end)


bins = np.linspace(start, end, int((end-start)/binning))

hist1, _ = np.histogram(time_stamps[0], bins)
hist2, _ = np.histogram(time_stamps[1], bins)
hist3, _ = np.histogram(time_stamps[2], bins)
hist4, _ = np.histogram(time_stamps[3], bins)

bins = bins[:-1]

print(len(bins), len(hist1))

fig, axs = plt.subplots(4, sharex = True, figsize = (6,10))

axs[0].plot(bins, hist1)
axs[0].set_title('a')
axs[1].plot(bins, hist2)
axs[1].set_title('b')
axs[2].plot(bins, hist3)
axs[2].set_title('c')
axs[3].plot(bins, hist4)
axs[3].set_title('d')

plt.title(f'Light Curves: {orbit}')

plt.savefig(f"{plot_path}/{orbit}/lc_{binning}.png")







