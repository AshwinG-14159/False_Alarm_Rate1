import csv
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle

import time
from scipy.ndimage import label

import subprocess
import glob

def fetch_files(source_path, destination_path, options=None):
    """
    Fetch files from source_path to destination_path using rsync.

    Parameters:
    - source_path: The source path for rsync.
    - destination_path: The destination path for rsync.
    - options: Additional options for rsync (optional).

    Returns:
    - Returns the return code of the rsync command.
    """
    rsync_command = ['rsync', '-av', '-e', 'ssh', source_path, destination_path]

    if options:
        rsync_command.extend(options)

    try:
        # Run rsync command
        result = subprocess.run(rsync_command, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        print(f"Error while running rsync: {e}")
        return e.returncode


def create_directory(directory_path):
    # Check if the directory exists
    if not os.path.exists(directory_path):
        # If it doesn't exist, create the directory
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")


t_script_start=time.time()

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


print("file reading time:",time.time()-t_script_start)
print("number of orbits to work with: ", len(set_of_rows))

path = "/mnt/nas2_czti/czti/"
path_tables = "../data/tables2"

t_comp_start = time.time()

Coinc_Tables = []
Analysed_Durations = []

skip_till = ['20230329_A12_054T02_9000005550_level2_40553', 'A12_054T02_9000005550', 'priyanka_iucaa', 'SBS 0846+513', '132.4916', '51.1414', '5100.90680462', '2023-03-29T16:16:53', '2023-03-29T18:16:54']
skipped = 1


err_counts = {"mkf file error":0, 
"coordinates error":0, 
"evt file header messup": 0, 
"evt file does not open":0,
"histograms not present":0,
"unable to make hists by region": 0}





orbit_num=0
for row in set_of_rows:
    orbit_num+=1
    if(row!=skip_till and not skipped):
        continue
    else:
        if(row==skip_till):
            skipped = 1
    t_row_start = time.time()
    print("orbit num: ", orbit_num)

    print('trying to access row:', row)
    date = row[0].split('_')[0]
    orbit = row[0].split('_')[-1]
    create_directory(f'../../data/evt_files/{orbit}')
    # target_folder = path+f'level2/{row[0][:-6]}/czti/orbit/{row[0][-5:]}_V1.0'
    # mkf_file = target_folder+f'/AS1{row[0][9:-13]}_{row[0][-5:]}czt_level2.mkf'
    # evt_file = target_folder+f'/modeM0/AS1{row[0][9:-13]}_{row[0][-5:]}cztM0_level2_bc.evt'
    
    # try:
    #     mkfdata = fits.getdata(mkf_file, 1)
    #     tab = Table.read(mkf_file,1)
    # except:
    #     err_counts["mkf file error"]+=1
    #     print("mkf file error")
    #     continue

    # Example usage
    files_needed = ['quad_clean.evt', 'quad_livetime.fits', 'quad_badpix.fits']
    source_path = f'cztipoc@192.168.11.37:/data2/czti/level2/{date}*level2/czti/orbit/{orbit}_V1.0/*czt_level2.mkf'#modeM0/*quad_badpix.fits'
    destination_path = f'../../data/evt_files/{orbit}/'

    return_code = fetch_files(source_path, destination_path, None)

    for file in files_needed:
        source_path = f'cztipoc@192.168.11.37:/data2/czti/level2/{date}*level2/czti/orbit/{orbit}_V1.0/modeM0/*{file}'
        destination_path = f'../../data/evt_files/{orbit}/'

        return_code = fetch_files(source_path, destination_path, None)


    # pattern_evt = f"*{orbit}*_clean.evt"


    # if return_code == 0:
    #     print("Files fetched successfully.")
    # else:
    #     print(f"Error fetching files. Return code: {return_code}")
    #     continue
    if(orbit_num==10):
        exit(0)
    else:
        continue

