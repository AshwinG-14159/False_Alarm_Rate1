import argparse
import glob

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import make_lightcurves_v2
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
destination_path = "/home/czti/user_area/ashwin/winter/data/evt_files"
pattern_quad_clean = f"*{orbit}*_quad_clean.evt"
pattern_mkf = f"*{orbit}*_level2.mkf"
pattern_badpix = f"*{orbit}*_quad_badpix.fits"
pattern_livetime = f"*{orbit}*_quad_livetime.fits"

plot_path = "../plots"
quad_clean_files = glob.glob(f'{destination_path}/*{pattern_quad_clean}')
mkf_files = glob.glob(f'{destination_path}/*{pattern_mkf}')
badpix_files = glob.glob(f'{destination_path}/*{pattern_badpix}')
livetime_files = glob.glob(f'{destination_path}/*{pattern_livetime}')

create_directory(f"{plot_path}/{orbit}")

time_stamps = []
quad_clean_file = quad_clean_files[0]
mkf_file = mkf_files[0]
badpix_file = badpix_files[0]
livetime_file = livetime_files[0]

for band in range(4):
    make_lightcurves_v2.bindata(quad_clean_file, mkf_file, badpix_file, livetime_file, 1, '', band)


exit(0)

orbitinfo, comb_startbins, comb_stopbins, comb_veto_startbins, comb_vetocift, comb_veto_stopbins, saa_start, saa_end = make_lightcurves_v2.make_detrended_lc(orbits, outpaths, args, dirname)





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
plt.cla()

l_curve, mask_lc,startbins3,stopbins3, error_flag = make_lightcurves_v2.getlc_clean(10, hist1, 5, 'median', 3, 5, 2)
print(l_curve.size, mask_lc.size, startbins3, stopbins3, error_flag)

plt.plot(bins,l_curve)
plt.savefig(f"{plot_path}/{orbit}/detrended_{binning}.png")
plt.cla()

