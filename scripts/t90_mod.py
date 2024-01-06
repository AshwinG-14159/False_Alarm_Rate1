from astropy.io import fits
import os.path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.stats import sigma_clipped_stats, sigma_clip
from matplotlib.backends.backend_pdf import PdfPages

import argparse

#------------------------------------------------------------------------
# Functions

def trend(time, a, b, c):
    """
    Given time series and some parameters, return the evaluated trend.
    """
    return a * time**2 + b * time + c


#------------------------------------------------------------------------
# Main code

def measure_t90(filename, plotfile, tmin, tmax, tran_start, tran_end, usequads, binsize, lc=[], tbins=[], newfilename='', newplotfile='', reject=True, lc_given=False, plotflag=True):
    if(lc_given == False):
        theplot = PdfPages(plotfile)
        quadnames = "".join([ ["A", "B", "C", "D"][quad] for quad in usequads ])
        suptitle = "Quad: {quadnames}, {binsize:0.2f}s binning, File:{filename}".format(quadnames=quadnames, binsize=binsize, filename=filename)
    else:
        theplot = newplotfile
        quadnames = "".join([ ["A", "B", "C", "D"][quad] for quad in usequads ])
        suptitle = "Quad: {quadnames}, {binsize:0.2f}s binning, File:{filename}".format(quadnames=quadnames, binsize=binsize, filename=newfilename)


    if(lc_given == False):
        d = None
        for quad in usequads:
            if d is None:
                d = fits.getdata(filename, quad+1)
            else:
                d = np.hstack((d, fits.getdata(filename, quad+1)))

        #d1 = fits.getdata(filename, 1)
        #d2 = fits.getdata(filename, 1)
        #d = np.hstack((d1, d2))
        #d = np.hstack((d1, d2, d3, d4))
        #d = np.copy(d2)

        d.sort(order='Time')
        select = (d['Time'] >= tmin) & (d['Time'] <= tmax)
        d = d[select]
        #print d, len(d)
        #tmin = min([min(d1['time']), min(d2['time'])])
        #tmax = max([max(d1['time']), max(d2['time'])])
        #tmin = np.min(d['Time'])
        #tmax = np.max(d['Time'])

        tbins = np.arange(tmin, tmax+binsize, binsize)
        lcurve, tbins = np.histogram(d['Time'], bins=tbins)
        lcurve = lcurve/binsize
    else:
        lcurve = lc
        tbins = tbins

    print(len(lcurve), len(tbins))
    if(len(tbins) == len(lcurve)):
        tmid = tbins
    else:
        tmid = 0.5 * (tbins[1:] + tbins[:-1])
    plotx = tmid - tmin
    mask = (tmid < tran_start) | (tmid > tran_end) # True for data points to include
    tranmask = np.copy(mask)

    if(plotflag):
        fig = plt.figure()
        plt.plot(tmid, lcurve, label="Data")
    #print 1
    p0 = [0, 0, np.median(lcurve[mask])]
    x = plotx[mask]
    popt, pcov = curve_fit(trend, x, lcurve[mask], p0)
    first_trend = trend(plotx, *popt)
    if(plotflag):
        plt.plot(tmid, first_trend, label="Initial trend")
    #print 2
    ratio = lcurve/first_trend # Exclude GRB while calculating ratio for sigma-clipping anyway
    rat_clipped = sigma_clip(ratio[mask], sigma=3)
    # rat_clipped.mask is False for data points to include in final fit
    clipmask = np.repeat(True, len(tmid))
    clipmask[mask] = ~rat_clipped.mask
    mask[mask] = ~rat_clipped.mask
    #print 3
    p0 = [0, 0, np.median(lcurve[mask])]
    x = plotx[mask]
    popt, pcov = curve_fit(trend, x, lcurve[mask], p0)
    final_trend = trend(plotx, *popt)
    if(plotflag):
        plt.plot(tmid, final_trend, label='Trend')
        plt.scatter(tmid[~tranmask], lcurve[~tranmask], marker='*', color='k', label='Transient')
        plt.scatter(tmid[~clipmask], lcurve[~clipmask], marker='o', color='r', label='Outliers')
        plt.legend(loc="upper left")
        plt.title("Actual data")
        plt.suptitle(suptitle)
        plt.axvspan(tran_start, tran_end, alpha=0.3, color='y')
        plt.xlim((tmin, tmax))
        plt.savefig('test_t90.png')
        theplot.savefig()
        plt.close(fig)
        fig = plt.figure()
        plt.plot(tmid, ratio)
        plt.title("Ratio")
        plt.suptitle(suptitle)
        plt.scatter(tmid[~tranmask], ratio[~tranmask], marker='*', color='k', label='Transient')
        plt.scatter(tmid[~clipmask], ratio[~clipmask], marker='o', color='r', label='Outliers')
        plt.axhline(1, linestyle='dashed', color='k')
        plt.legend(loc="upper left")
        plt.axvspan(tran_start, tran_end, alpha=0.3, color='y')
        plt.xlim((tmin, tmax))
        theplot.savefig()
        plt.close(fig)
    #print 4
    
    detrend_lc = lcurve - final_trend
    if(plotflag):
        fig = plt.figure()
        plt.plot(tmid, detrend_lc)
        plt.title("Detrended LC")
        plt.suptitle(suptitle)
        plt.axvspan(tran_start, tran_end, alpha=0.3, color='y')
        plt.axhline(0, linestyle='dashed', color='k')
        plt.xlim((tmin, tmax))
        theplot.savefig()
        plt.close(fig)
        fig = plt.figure()
        plt.title(r"T$_{90}$ measurement")
        plt.suptitle(suptitle)
    cs0 = np.cumsum(detrend_lc) / np.sum(detrend_lc)
    cs = np.cumsum(detrend_lc[clipmask]) / np.sum(detrend_lc[clipmask])
    #print 5
    ## get background and rate above background
    bg = np.mean(final_trend[~mask])
    rate = np.max(detrend_lc[~mask])
    counts = np.sum(detrend_lc[~mask])*binsize
    #print(bg,rate,counts)
    #print 6
    if reject:
        # Reject outliers in T90 calculation
        #print "Rejecting sigma clip outliers"
        #print "Total number of photons: ", np.sum(detrend_lc[clipmask])
        tclip = tmid[clipmask]
        #selclip = (tclip >= tran_start) & (tclip <= tran_end)
        selclip = tclip > 0
        t_trans = tclip[selclip]
        cs_trans = cs[selclip]
    else:
        # No outlier rejection
        #print "Keeping sigma clip outliers"
        t_trans = tmid[~tranmask]
        cs_trans = cs[~tranmask]

    #print 7    
    # Measure T90 in the transient region
    #cs_check = np.array([0.05, 0.95])
    #t_check = np.interp(cs_check, cs_trans, t_trans) ## Interpolation fails with mutliple intersections
    index_05 = np.max(np.where(cs_trans <= 0.05))
    t05 = np.interp(0.05, cs_trans[index_05:index_05+2], t_trans[index_05:index_05+2]) # Last element of range is not included
    index_95 = np.min(np.where(cs_trans >= 0.95))
    if(index_95 == 0):
        index_95 = 1
    t95 = np.interp(0.95, cs_trans[index_95-1:index_95+1], t_trans[index_95-1:index_95+1]) # Last element of range is not included
    t90 = t95 - t05
    print("T90 is {t90:0.1f} for binsize {binsize:0.2f}".format(t90=t90, binsize=binsize))
    t90text = r"T$_{{90}}$ = {t90:0.2f}".format(t90=t90)

    if(plotflag):
        plt.plot(tmid, cs0, label="Without outlier rejection")
        plt.plot(tmid[clipmask], cs, label="With outlier rejection")
        plt.axhline(0.05, color='gray', alpha=0.5)
        plt.axhline(0.95, color='gray', alpha=0.5)
        plt.axhline(0.00, color='k')
        plt.axhline(1.00, color='k')
        plt.axvline(t05, color='r', label=t90text)
        plt.axvline(t95, color='r')
        #plt.scatter(t_trans, cs_trans)
        plt.axvspan(tran_start, tran_end, alpha=0.3, color='y')
        plt.legend(loc="best")
        plt.xlim((tmin, tmax))
        theplot.savefig()
        plt.close(fig)

    if(lc_given == False):
        theplot.close()
    #print "Plots saved to ", plotfile
    print(bg, rate, counts)
    return t90, bg, rate, counts

def iterate_t90(filename, plotfile, tmin, tmax, tran_start, tran_end, usequads, binsize, lc, tbins, newfilename, newplotfile, reject=True, lc_given=True):
    prev_t90 = 0
    count = 0
    while True:
        print(tran_start, tran_end)
        try:
            if(count==5):
                plotflag=True
            else:
                plotflag=False
            t90, bg, rate, counts = measure_t90(filename, plotfile, tmin, tmax, tran_start, tran_end, usequads, binsize, lc, tbins, newfilename, newplotfile, reject, lc_given, plotflag)
        except:
            t90, bg, rate, counts = -1, 0, 0, 0
        tran_start = tran_start - binsize
        tran_end = tran_end + binsize        
        if(count == 5):
            break
        prev_t90 = t90
        count = count + 1
    #print t90, bg, rate, counts
    return t90, bg, rate, counts




#------------------------------------------------------------------------

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(epilog="""
        Measure the T90 of a GRB by using various bin sizes.
        The duration from transtart to tranend must definitely be completely
        contained with the GRB: which also means T90 calculated by this code
        will always be longer than 'tranend-transtart'.
        """)
    parser.add_argument("evtfile", help="An event file to be processed", type=str)
    parser.add_argument("plotbase", nargs='?', type=str, help="Stem for plots showing noise-dominated times. Don't include an extension", default="products/t90plots")
    parser.add_argument("--tmin", help="Start time (in spacecraft seconds). Default = start of file", type=float, default=None)
    parser.add_argument("--tmax", help="End time (in spacecraft seconds). Default = end of file", type=float, default=None)
    parser.add_argument("--transtart", nargs=1, help="First instant where GRB emission is very clearly seen (in spacecraft seconds).", type=float, default=None)
    parser.add_argument("--tranend", nargs=1, help="Last instant where GRB emission is very clearly seen (in spacecraft seconds).", type=float, default=None)
    parser.add_argument("--quad", type=str, help="Quadrants to be processed (for instance ABD, default=ABCD", default="ABCD")
    parser.add_argument("--binsize", nargs='*', help="Bin sizes to be used while calculating T90",
                        default=[0.1, 0.3, 1.0, 3.0],type=float)
    parser.add_argument("--noreject", dest='reject', action='store_false') # If option not given, args.reject becomes True
    parser.add_argument("-q", "--quiet", dest='verbose', action='store_false') # If option not given, args.verbose becomes True
    args = parser.parse_args()
    #print args
    
    
    hdu = fits.open(args.evtfile)
    # If needed, calculate tmin and tmax
    if args.tmin is None:
        #if args.verbose: print "Determining tmin and tmax"
        args.tmin = 1e15
        for quad in range(4):
            if len(hdu[quad+1].data) > 0:
                newmin = min(hdu[quad+1].data['Time'])
                args.tmin = min([args.tmin, newmin])

    if args.tmax is None:
        #if args.verbose: print "Determining tmin and tmax"
        args.tmax = 0.
        for quad in range(4):
            if len(hdu[quad+1].data) > 0:
                newmax = max(hdu[quad+1].data['Time'])
                args.tmax = max([args.tmax, newmax])

    if (args.transtart is None) or (args.tranend is None):
        print("Error! Both 'transtart' and 'tranend' must be specified!")
        raise SystemExit

    # Figure out which quadrants to use
    usequads = []
    quadlist = list(args.quad.upper())

    for quadnum, quadname in enumerate(["A", "B", "C", "D"]):
        if quadname in quadlist: usequads.append(quadnum)

    #bisizes =[0.07, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00, 1.50, 2.00, 2.50, 3.00]
    for binsize in args.binsize:
        plotfile = "{stem}{bins:0.3f}.pdf".format(stem=args.plotbase, bins=binsize)
        measure_t90(args.evtfile, plotfile, args.tmin, args.tmax, args.transtart[0], args.tranend[0], usequads, binsize, reject=args.reject,plotflag=True)
