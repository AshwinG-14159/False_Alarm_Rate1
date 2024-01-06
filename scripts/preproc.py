import os
import copy
import subprocess
import argparse
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.signal import savgol_filter
import datetime
from astropy.table import Table, Column, vstack
# import astrosat_time
# import plot
import pandas as pd

def gethist(lc_all, masks_all):
	"""Get histogram of rates for all quadrants.

	Inputs:
	lc_all - List of light curves of all quadrants
	masks_all - List of all masks

	Outputs:
	hist_all - List of [rate_bins, rate_hist] for each quadrants

	Local variables:
	ratesbin - Bin size of rates
	min_rate - minimum rate for each quadrants
	max_rate - maximum rate for each quadrant
	rate_hist, rate_bins - histogram and bins of rates
	"""
	ratesbin = 0.1
	hist_all = []
	#print masks_all
	for quadno in range(4):           
		if(len(lc_all[quadno][~masks_all[quadno]])==0):
			fakelc=np.zeros((100))
			rate_hist, rate_bins = np.histogram(fakelc, bins=np.arange(0,10+0.1,0.1))
			hist_all.append([rate_bins,rate_hist])
		else:
			min_rate = np.min(lc_all[quadno][~masks_all[quadno]])
			max_rate = np.max(lc_all[quadno][~masks_all[quadno]])
			if(min_rate < -1000):
				min_rate = -1000
			bins = np.arange(min_rate, max_rate + ratesbin, ratesbin)
			rate_hist, rate_bins = np.histogram(lc_all[quadno][~masks_all[quadno]],bins=bins,range=(bins.min(),bins.max()))
			hist_all.append([rate_bins, rate_hist])
	return hist_all

def getcutoffrate1(tbin, hist_all, far, timespan,lc_all,mask_all):
	"""Get cutoff rate for each quadrant for a given quadrant.

	False Alarm Rate 'far' corresponds to a confidence of
	1.0 - tbin*far/timespan
	Also, confidence = cummulative sum/total sum. The cutoffrate can be found
	by equating the two equations

	Inputs:
	tbin - bin_times
	hist_all - Histogram of rates for each quadrant - [rate_bins, rate_hist]
	far - False alarm rate
	timespan - Time over which the far is to be evaluated
	stats - Gaussian distribution parameters from sigma clip stats

	Outputs:
	cutoffs - List of cutoffrates for each quadrant

	Local variables:
	confidence - confidence corresponding to given far
	quadno - quadrant number
	rate_bins - Mid value of the time bins
	cum_frac - Cummulative sum fraction
	cut_index, cut_findex - indices corresponding to the confidence in cum_frac
	cut_rate - cutoff rate for the given quadrant
	"""
	confidence = 1.0 - tbin*far/timespan
	cutoffs = []
	for quadno in range(4):
		#calculating stats
		#stats = sigma_clipped_stats(lc_all[quadno],mask=mask_all[quadno],sigma=5.0,iters=10)
		#gaussian = np.exp(-0.5*np.power(((hist_all[quadno][0]-stats[0])/stats[2]),2.))/(2*np.pi*(stats[2]**2))
		# Taking cummulative sum
		cum_frac = (np.cumsum(hist_all[quadno][1])*1.0 /
					np.sum(hist_all[quadno][1]))
		if(len(cum_frac)==0):
			cutoffs.append(0)
			continue
		#cum_frac = (np.cumsum(gaussian)*1.0 /
		#            np.sum(gaussian))
		# Searching for the index corresponding to confidence
		cut_index = np.searchsorted(cum_frac, confidence)
		if(cut_index == 0):
			cut_index = 1
		# Getting the fractional part of the index where cummulative sum is
		# equal to confidence
		cut_findex = np.interp(confidence, cum_frac[cut_index-1:cut_index+1],[cut_index, cut_index+1])
		cut_rate = np.interp(cut_findex, [cut_index, cut_index+1],hist_all[quadno][0][cut_index:cut_index+2])
		cutoffs.append(cut_rate)
	return cutoffs

def getcutoffrate2(lc_all, mask_all, nsigma):
	"""Get cutoff rates based on nsigma cutoff.

	Clip the data using astropy sigma clipping.

	Inputs:
	lc_all - Light curves of each quadrant
	mask_all - masks in each light curve
	nsigma - cutoff number of sigmas

	Outputs:
	cutoffs - cut off rates in each quadrant
	"""
	cutoffs = []
	for quadno in range(4):
		stats = sigma_clipped_stats(lc_all[quadno], mask=mask_all[quadno],
									sigma=5.0, maxiters=3)
		cutrate = stats[0] + nsigma*stats[2]
		cutoffs.append(cutrate)
	return cutoffs


def getpeaks(lc_all, cutoffs):
	"""Return times when peak has occured in atleast 2 quadrants.

	Gets a binary map with 1 where rate > cutoffrate for each quadrant.
	Add the 4 maps and then get the times where sum > 2

	Inputs:
	lc_all - light curves of all quadrants
	cutoffs - cutoff rates in each quadrant

	Returns
	peak_map - Binary map of the peak bins
	positions where the peaks are coinicident in atleast 2 quadrants

	"""
	arr_len = max(len(lc_all[0]),len(lc_all[1]),len(lc_all[2]),len(lc_all[3]))
	peak_map = np.zeros((4,arr_len))
	lc_all = np.array(lc_all)
	cutoffs = np.array(cutoffs)
	#print len(lc_all[0]),len(lc_all[1]),len(lc_all[2]),len(lc_all[3])
	peak_map[lc_all > cutoffs[:,None]] = 1
	# for quad in range(4):
	# 	for index in range(arr_len):
	# 		lc_element = lc_all[quad][index]
	# 		cutoff = cutoffs[quad]
	# 		#print lc_element, cutoff
	# 		if(lc_element > cutoff):
	# 			peak_map[quad][index] = 1
	#print peak_map
	return np.where(np.sum(peak_map, axis=0) >= 2.0), peak_map

def getpeakstopN(lc_all,N):
	"""Return times when peak has occured in atleast 2 quadrants.

	Gets a binary map with 1 where rate > cutoffrate for each quadrant.
	Add the 4 maps and then get the times where sum > 2

	Inputs:
	lc_all - light curves of all quadrants
	cutoffs - cutoff rates in each quadrant

	Returns
	peak_map - Binary map of the peak bins
	positions where the peaks are coinicident in atleast 2 quadrants (not using from here)

	"""
	mask = np.zeros((4,len(lc_all[0])))
	lc_all = np.array(lc_all)
	inds = lc_all.argsort()[:,-3:][:,::-1]
	for q in range(4):
		mask[q][inds[q]] = 1
	peak_map = mask
	return np.where(np.sum(peak_map, axis=0) >= 2.0), peak_map 


# def getconf(lc_all, peak_inds, peak_map, mask_all):
# 	"""Get total confidence in the detected GRBs.

# 	Inputs:
# 	lc_all - list of all light curves
# 	peak_ind - Indices where the light curve peaks
# 	peak_map - Binary map of peaks
# 	mask_all - refined mask after removing peaks

# 	Outputs:
# 	conf - confidence in each of the quadrants
# 	conf_al - total confidence across 4 quadrants
# 	conf_det - total confience in the detected quadrants
# 	"""
# 	#print peak_inds
# 	conf = np.zeros((len(peak_inds), 4))
# 	conf_det = np.zeros(len(peak_inds))
# 	for peakno, peak in enumerate(peak_inds):
# 		for quadno in range(4):
# 			conf[peakno, quadno] = (
# 				len(np.where(lc_all[quadno][~mask_all[quadno]] <
# 							 lc_all[quadno][peak])[0])*1.0 /
# 				len(lc_all[quadno][~mask_all[quadno]]))
# 		det_quad = np.where(peak_map[:, peak] == 1)[0]
# 		conf_det[peakno] = np.prod(conf[peakno, det_quad])
# 	conf_all = np.prod(conf, axis=1)
# 	return conf, conf_all, conf_det

def grb_lc(combined_lc, combined_mask, combined_tbins, grb_peakind, tbin):
	num_bins = 50 
	num_bins = int(num_bins)
	minbin = max(0,grb_peakind-num_bins)
	maxbin = min(len(combined_tbins),grb_peakind+num_bins+1)
	#print len(combined_tbins),len(combined_lc[0]),minbin,maxbin
	tplot = combined_tbins[minbin:maxbin]
	combined_lc = np.array(combined_lc)
	lc = combined_lc[:,:,minbin:maxbin-1]
	mask=np.zeros(np.shape(lc))
	for p in range(3):
		for q in range(4):
			mask[p][q]=combined_mask[p][q][minbin:maxbin-1]
	# print"-------------"
	# print np.shape(mask), np.shape(lc)
	#print len(tplot),len(lc[0][0])
	return tplot,lc,mask

def total_countrate(bins_time,peakinds,peak_ind_bin,combined_lc):
	total_counts=[]
	#print peakinds
	for peakind in peakinds:
		count=combined_lc[:,:,peakind]
		# for x in range(len(peak_ind_bin)):
		# 	if(bins_time[peak_ind_bin[x]]-bins_time[peakind]<100.0 and bins_time[peak_ind_bin[x]]-bins_time[peakind]>=0):
		# 		count=count+combined_lc[:,:,peak_ind_bin[x]]
		total_counts.append(count)
	return total_counts

def get_rank(quads):
	rank=0
	band0=[]
	band1=[]
	band2=[]
	for num, quad in enumerate(quads):
		if(quad>= 0 and quad<=3):
			band0.append(quad)
			rank=rank+1
		elif(quad>=4 and quad<=7):
			band1.append(quad-4)
			rank=rank+2
		elif(quad>=8 and quad<12):
			band2.append(quad-8)
			rank=rank+3

	return rank,band0,band1,band2

def get_significance(lc,mask,total_counts):
	stat_significance = []
	band_significance = []
	for band in range(3):
		for quad in range(4):
			stats = sigma_clipped_stats(lc[band][quad],mask=mask[band][quad], sigma = 3.0, maxiters = 2)
			if(stats[2]!=0):
				sig = (max(lc[band][quad])-stats[0])/stats[2]
			else:
				sig = (max(lc[band][quad])-stats[0])
			stat_significance.append(sig)
		tot_sig = (stat_significance[4*band+0]**2 + stat_significance[4*band+1]**2 + stat_significance[4*band+2]**2 + 
					stat_significance[4*band+3]**2)**0.5
		band_significance.append(tot_sig)
	return stat_significance,band_significance

def total_countrate_veto(bins_time,peakinds,peak_ind_bin,combined_lc):
	total_counts=[]
	for peakind in peakinds:
		count=combined_lc[:,peakind]
		# for x in range(len(peak_ind_bin)):
		# 	if(bins_time[peak_ind_bin[x]]-bins_time[peakind]<100.0 and bins_time[peak_ind_bin[x]]-bins_time[peakind]>=0):
		# 		count=count+combined_lc[:,peak_ind_bin[x]]
		total_counts.append(count)
	return total_counts

def grb_lc_veto(combined_lc,combined_mask, combined_tbins,grb_peakind,tbin):
	num_bins = 30 
	num_bins = int(num_bins)
	minbin = max(0,grb_peakind-num_bins)
	maxbin = min(len(combined_tbins),grb_peakind+num_bins+1)
	#print len(combined_tbins),len(combined_lc[0]),minbin,maxbin
	tplot = combined_tbins[minbin:maxbin]
	lc = combined_lc[:,minbin:maxbin]
	mask=np.zeros(np.shape(lc))
	for q in range(4):
		mask[q]=combined_mask[q][minbin:maxbin]	
	#print len(tplot),len(lc[0][0])
	return tplot,lc,mask

def get_significance_veto(lc,mask,total_counts):
	stat_significance = []
	for quad in range(4):
		stats = sigma_clipped_stats(lc[quad],mask=mask[quad], sigma = 3.0, maxiters = 2)
		if(stats[2]!=0):
			sig = (max(lc[quad])-stats[0])/stats[2]
		else:
			sig = max(lc[quad])-stats[0]
		stat_significance.append(sig)
	if(np.any(stat_significance) == np.inf or np.any(np.isnan(stat_significance))):
		tot_sig = 0.0
	else:
		tot_sig = (stat_significance[0]**2 + stat_significance[1]**2 + stat_significance[2]**2 + 
						stat_significance[3]**2)**0.5
	return stat_significance,tot_sig

def plot_hist(hists,cutoffs,title,pdffile):
	for quad in range(4):
		bins = hists[quad][0]
		bins = 0.5*(bins[1:]+bins[:-1])
		plt.plot(bins,hists[quad][1],label='Quad '+str(quad))
		plt.show()
		plt.axvline(cutoffs[quad],color='black',linestyle='dashed')
		plt.title(title)
		plt.xlabel('Rates')
		plt.ylabel('Rate histogram')
		pdffile.savefig()
		plt.close()