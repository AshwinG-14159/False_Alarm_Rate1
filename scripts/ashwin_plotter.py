import argparse
import glob

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import make_lightcurves_v2
import os
import algorithms
from astropy.table import Table, Column, vstack
from astropy.stats import sigma_clipped_stats
import csv

# def makeparser():
#         """Make parser."""
#         parser = argparse.ArgumentParser(
#                 description="""Program to process AstroSAT data and detect GRBs.
#                 Uses the czti pipeline to process the data and get light curves.
#                 Applies median filter to filter the data.
#                 Compare the peaks across quadrants and detect GRBs
#                 """,
#                 epilog="""Cutoff rate calculated from false alam rates.
#                 Last Update:
#                 """)
#         """parser.add_argument('infolder', help='Folder(orbit number) on which' +
#                                                 'analysis needs to be done', type=str)
#         parser.add_argument('outpath', nargs='?', default=None,
#                                                 help='Path for the creation of output' +
#                                                 'folders and files', type=str)
#         parser.add_argument("outbase", nargs="?", default=None, type=str,
#                                                 help="Stem to be prepended to all output files.") """
#         parser.add_argument('directories', help="Runs on which analysis needs to. Select from the given list. Default is none"
#                                                 , nargs='+',type=str)
#         parser.add_argument('-p','--pipeline', type=str, choices=['yes','no','auto'],default='auto',
#                                                 help='Select "yes" to run pipeline, "no" to skip running pipeline'+
#                                                 'and "auto" if confused or do not care. By default I will assume you are confused :P')
#         # parser.add_argument('-g','--grb',type=str,choices=['yes','no'],default='yes',
#         # 					help='Select "yes" to run grb search, "no" to skip')
#         parser.add_argument('-f', '--far', help='False alarm rate per quadrant' +
#                                                 'Default = 1', type=float, default=1.0)
#         parser.add_argument("--timespan", "--ts", "--roi", help="Time span of" +
#                                                 "transient search. This is the region of interest" +
#                                                 "used in estimating FAR, default=100",
#                                                 type=float, default=4156)
#         parser.add_argument("--filtertype", help="'median' or 'savgol'", type=str,
#                                                 default='savgol')
#         parser.add_argument("-w", "--filterwidth", help="Filterwidth in seconds" +
#                                                 "Default=100", type=float, default=100)
#         parser.add_argument("-o", "--filter_order", "--order", help="Polynomial" +
#                                                 "order for SG filter, default=2", type=int, default=2)
#         parser.add_argument("-b", "--tbin", help="Bin sizes for output lc." +
#                                                 "Default = [0.1, 1.0, 10.0]", action="append",
#                                                 type=float)
#         parser.add_argument("-t", "--threshold", help="Threshold count rate for" +
#                                                 "finding contiguous windows", type=float, default=0.0)
#         parser.add_argument("--tclip", default=15.0, type=float, help="Time in" +
#                                                 "seconds to be clipped at the begining and at the" +
#                                                 "end of window after detrending")
#         # parser.add_argument("--cutoffmethod",type=str,choices=['nsigma','cumsum','topN'],default='cumsum',
#         # 					help='Select the type of method for getting cutoff rates for peak selection')
#         parser.add_argument("-plot", "--plottype", choices=['png','pdf'], help="Type of output file to produce (png, pdf)", type=str, default='pdf')
#         # parser.add_argument("--runtype",type=str,choices=['orbit','obs-id'],default='obs-id',
#         # 					help="Option to run cutoff algorithm on each orbit or whole obs-id")
#         parser.add_argument('-pf','--plotflag', action='store_true')
#         parser.add_argument('-mvrm','--moveandremove',action='store_false',help='specify this argument to disable moving and removing of data files')
#         return parser



def create_directory(directory_path):
    # Check if the directory exists
    if not os.path.exists(directory_path):
        # If it doesn't exist, create the directory
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")


def plot_quads(times_arr, rates_arr, heading, filename):
    fig, axs = plt.subplots(4, sharex = True, figsize = (6,10))

    axs[0].plot(times_arr[0], rates_arr[0])
    axs[0].set_title('a')
    axs[1].plot(times_arr[1], rates_arr[1])
    axs[1].set_title('b')
    axs[2].plot(times_arr[2], rates_arr[2])
    axs[2].set_title('c')
    axs[3].plot(times_arr[3], rates_arr[3])
    axs[3].set_title('d')

    plt.suptitle(f'{heading}')

    plt.savefig(f"{plot_path}/{orbit}/{filename}.png")
    plt.cla()
    plt.close()

def give_window_lcs(startbins, stopbins, times, rates):
    times_arr = []
    rates_arr = []
    for i in range(len(startbins)):
        # print('s',startbins[i], 's',stopbins[i])
        cond = (times>=startbins[i])&(times<=stopbins[i])
        # print(cond)
        times_arr.append(times[cond])
        rates_arr.append(rates[cond])
    return times_arr, rates_arr



def examine_stats(rates, times):

    pass

set_of_rows = []

gap = 1
total_attempts = 100

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


    binnings = [0.1,1,10]

    my_path = f"../../data/evt_files/{orbit}"


    plot_path = "../plots"
    result_path = '../results'

    log_path = "../logs"



    # create_directory(f"{plot_path}/{orbit}")
    # create_directory(f"{log_path}/{orbit}")
    # create_directory(f"{result_path}/{orbit}")


    # pattern_lc = f"*.lc"
    # lc_files = glob.glob(f'{my_path}/*{pattern_lc}')
    # print(lc_files)
    # parser = makeparser()
    # args = parser.parse_args()
    # args.tbin = binnings



    # orbitinfo, comb_startbins, comb_stopbins, comb_veto_startbins, comb_veto_stopbins, saa_start, saa_end = make_lightcurves_v2.make_detrended_lc(orbits = [my_path], outpaths=[f'{my_path}'], args=args, dirname='.')


    threshs = []
    counts = []

    for binning_id in range(len(binnings)):
        binning = binnings[binning_id]
        result_file = np.load(f'{result_path}/{orbit}/{orbit}_{binning}_nsigma_FAR.npy')
        thresh_arr = result_file[0,:][1:]
        count_arr = result_file[1,:][1:]/result_file[1,1]
        threshs.append(thresh_arr)
        count.append(count_arr)


        # print(list_of_counts)
    # for i in range
        plt.plot(thresh_arr, count_arr, drawstyle = "steps", label = f'{binning}')
    
    plt.title(f"False Alarm Rates, Orbit {orbit}, binning: {binning}, alg:nsigma")
    plt.yscale('log')
    plt.xlabel("nsigma threshold")
    plt.ylabel("False Alarm Rate")
    plt.savefig(f"{plot_path}/{orbit}/{orbit}_False_alarm_counts_combined.png")
    plt.close()




    # exit(0)

