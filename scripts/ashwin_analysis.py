import argparse
import glob

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import make_lightcurves_v2
import os
import algorithms

def makeparser():
        """Make parser."""
        parser = argparse.ArgumentParser(
                description="""Program to process AstroSAT data and detect GRBs.
                Uses the czti pipeline to process the data and get light curves.
                Applies median filter to filter the data.
                Compare the peaks across quadrants and detect GRBs
                """,
                epilog="""Cutoff rate calculated from false alam rates.
                Last Update:
                """)
        """parser.add_argument('infolder', help='Folder(orbit number) on which' +
                                                'analysis needs to be done', type=str)
        parser.add_argument('outpath', nargs='?', default=None,
                                                help='Path for the creation of output' +
                                                'folders and files', type=str)
        parser.add_argument("outbase", nargs="?", default=None, type=str,
                                                help="Stem to be prepended to all output files.") """
        parser.add_argument('directories', help="Runs on which analysis needs to. Select from the given list. Default is none"
                                                , nargs='+',type=str)
        parser.add_argument('-p','--pipeline', type=str, choices=['yes','no','auto'],default='auto',
                                                help='Select "yes" to run pipeline, "no" to skip running pipeline'+
                                                'and "auto" if confused or do not care. By default I will assume you are confused :P')
        # parser.add_argument('-g','--grb',type=str,choices=['yes','no'],default='yes',
        # 					help='Select "yes" to run grb search, "no" to skip')
        parser.add_argument('-f', '--far', help='False alarm rate per quadrant' +
                                                'Default = 1', type=float, default=1.0)
        parser.add_argument("--timespan", "--ts", "--roi", help="Time span of" +
                                                "transient search. This is the region of interest" +
                                                "used in estimating FAR, default=100",
                                                type=float, default=4156)
        parser.add_argument("--filtertype", help="'median' or 'savgol'", type=str,
                                                default='savgol')
        parser.add_argument("-w", "--filterwidth", help="Filterwidth in seconds" +
                                                "Default=100", type=float, default=100)
        parser.add_argument("-o", "--filter_order", "--order", help="Polynomial" +
                                                "order for SG filter, default=2", type=int, default=2)
        parser.add_argument("-b", "--tbin", help="Bin sizes for output lc." +
                                                "Default = [0.1, 1.0, 10.0]", action="append",
                                                type=float)
        parser.add_argument("-t", "--threshold", help="Threshold count rate for" +
                                                "finding contiguous windows", type=float, default=0.0)
        parser.add_argument("--tclip", default=15.0, type=float, help="Time in" +
                                                "seconds to be clipped at the begining and at the" +
                                                "end of window after detrending")
        # parser.add_argument("--cutoffmethod",type=str,choices=['nsigma','cumsum','topN'],default='cumsum',
        # 					help='Select the type of method for getting cutoff rates for peak selection')
        parser.add_argument("-plot", "--plottype", choices=['png','pdf'], help="Type of output file to produce (png, pdf)", type=str, default='pdf')
        # parser.add_argument("--runtype",type=str,choices=['orbit','obs-id'],default='obs-id',
        # 					help="Option to run cutoff algorithm on each orbit or whole obs-id")
        parser.add_argument('-pf','--plotflag', action='store_true')
        parser.add_argument('-mvrm','--moveandremove',action='store_false',help='specify this argument to disable moving and removing of data files')
        return parser



def create_directory(directory_path):
    # Check if the directory exists
    if not os.path.exists(directory_path):
        # If it doesn't exist, create the directory
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")


def plot_quads(times_arr, rates_arr, orbit, filename):
    fig, axs = plt.subplots(4, sharex = True, figsize = (6,10))

    axs[0].plot(times_arr[0], rates_arr[0])
    axs[0].set_title('a')
    axs[1].plot(times_arr[1], rates_arr[1])
    axs[1].set_title('b')
    axs[2].plot(times_arr[2], rates_arr[2])
    axs[2].set_title('c')
    axs[3].plot(times_arr[3], rates_arr[3])
    axs[3].set_title('d')

    plt.suptitle(f'Light Curves: {orbit}')

    plt.savefig(f"{plot_path}/{orbit}/{filename}.png")
    plt.cla()




def examine_stats(quad_clean):
    pass



orbit = 43416
binning = 10
destination_path = "/home/czti/user_area/ashwin/winter/data/evt_files"
my_path = "../../data/evt_files"
output_path = "/home/czti/user_area/ashwin/winter/data/outputs"

# pattern_quad_clean = f"*{orbit}*_quad_clean.evt"
# pattern_mkf = f"*{orbit}*_level2.mkf"
# pattern_badpix = f"*{orbit}*_quad_badpix.fits"
# pattern_livetime = f"*{orbit}*_quad_livetime.fits"


plot_path = "../plots"




# quad_clean_files = glob.glob(f'{destination_path}/*{pattern_quad_clean}')
# mkf_files = glob.glob(f'{destination_path}/*{pattern_mkf}')
# badpix_files = glob.glob(f'{destination_path}/*{pattern_badpix}')
# livetime_files = glob.glob(f'{destination_path}/*{pattern_livetime}')

# create_directory(f"{plot_path}/{orbit}")

# time_stamps = []
# quad_clean_file = quad_clean_files[0]
# mkf_file = mkf_files[0]
# badpix_file = badpix_files[0]
# livetime_file = livetime_files[0]

pattern_lc = f"*.lc"
lc_files = glob.glob(f'{my_path}/*{pattern_lc}')
# print(lc_files)
parser = makeparser()
args = parser.parse_args()
args.tbin = [1]



orbitinfo, comb_startbins, comb_stopbins, comb_veto_startbins, comb_veto_stopbins, saa_start, saa_end = make_lightcurves_v2.make_detrended_lc(orbits = [my_path], outpaths=[my_path], args=args, dirname='.')

for band in range(3):
    quad_lcs_rates = []
    quad_lcs_times = []
    for quad in range(4):
        lc_file = f"{destination_path}/_1_{band}_Q{quad}.lc"
        with fits.open(lc_file) as hdul:
            quad_lcs_rates.append(hdul[1].data['rate'])
            quad_lcs_times.append(hdul[1].data['time'])
    
    plot_quads(quad_lcs_times, quad_lcs_rates, orbit, f"lc_{binning}_{band}")

for band in range(3):
    quad_lcs_rates = []
    quad_lcs_times = []
    for quad in range(4):
        lc_file = f"{destination_path}/_1_{band}_Q{quad}.lc"
        with fits.open(lc_file) as hdul:
            quad_lcs_rates.append(hdul[1].data['rate'])
            quad_lcs_times.append(hdul[1].data['time'])
    
    plot_quads(quad_lcs_times, quad_lcs_rates, orbit, f"lc_{binning}_{band}")


numgrb_nsigma = algorithms.grb_search(args, [my_path], orbitinfo, comb_startbins, comb_stopbins, '.', 'nsigma', saa_start, saa_end)

print(numgrb_nsigma)












exit(0)






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

