import argparse
import glob
import csv

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


set_of_rows = []

binning = 0.1
gap = 1
total_attempts = 100
version = 1

with open('../../data/orbitinfo.csv', 'r') as f:
    r = csv.reader(f)
    count=0
    for row in r:

        if(count<1):
            count+=1
            continue

        if(float(row[5])>-180 and float(row[5])<180) and row[0][-6:]!='level2':
            if(count%gap!=0):
                count+=1
                continue
            set_of_rows.append(row)
            count+=1
            if(count/gap>total_attempts): # leave 10 between any 2 and take a total of 200. Binning = 0.001
                break

skip_till = ['20230329_A12_054T02_9000005550_level2_40553', 'A12_054T02_9000005550', 'priyanka_iucaa', 'SBS 0846+513', '132.4916', '51.1414', '5100.90680462', '2023-03-29T16:16:53', '2023-03-29T18:16:54']
skipped = 1
orbit_num=0
for row in set_of_rows:
    orbit_num+=1
    if(row!=skip_till and not skipped):
        continue
    else:
        if(row==skip_till):
            skipped = 1
    # t_row_start = time.time()
    print("orbit num: ", orbit_num)

    print('trying to access row:', row)
    date = row[0].split('_')[0]
    orbit = row[0].split('_')[-1]
    date = "20190926"
    orbit = "21599"


    # orbit = 43416
    binning = [0.064, 0.032, 0.016, 0.1, 1, 10]
    destination_path = f"/home/czti/user_area/ashwin/winter/data/evt_files/{orbit}/"
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
    # print(quad_clean_files, f'{destination_path}/*{pattern_quad_clean}')
    try:
        quad_clean_file = quad_clean_files[0]
        mkf_file = mkf_files[0]
        badpix_file = badpix_files[0]
        livetime_file = livetime_files[0]
    except:
        continue

    for bin in binning:
        for band in range(3):
            make_lightcurves_v2.bindata(quad_clean_file, mkf_file, badpix_file, livetime_file, bin, f'{orbit}', band)
            # exit(0)
    exit(0)
    if(orbit_num==120):
        exit(0)
    else:
        continue

