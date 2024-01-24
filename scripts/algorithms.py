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
import preproc
import sqlite3

def add_czti_event(mdb, event):
	sql = ''' INSERT OR IGNORE INTO events (eventid, triggertime, binsize, method, orbit, obsid, preSAA, postSAA, rank, quadsband0, quadsband1, 
	quadsband2, significanceband0, significanceband1, significanceband2, sigb0Q0, sigb0Q1, sigb0Q2, sigb0Q3,
	sigb1Q0, sigb1Q1, sigb1Q2, sigb1Q3, sigb2Q0, sigb2Q1, sigb2Q2, sigb2Q3, rateb0Q0, rateb0Q1, rateb0Q2, 
	rateb0Q3, rateb1Q0, rateb1Q1, rateb1Q2, rateb1Q3, rateb2Q0, rateb2Q1, rateb2Q2, rateb2Q3) 
	VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) '''
	cur = mdb.cursor()
	print(event)
	cur.execute(sql, event)
	return cur.lastrowid

def add_veto_event(mdb, event):
	sql = ''' INSERT OR IGNORE INTO events (eventid, triggertime, binsize, method, orbit, obsid, preSAA, postSAA, quadsband0, significanceband0, 
	sigb0Q0, sigb0Q1, sigb0Q2, sigb0Q3, rateb0Q0, rateb0Q1, rateb0Q2, rateb0Q3) 
	VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) '''
	cur = mdb.cursor()
	print(event)
	cur.execute(sql, event)
	return cur.lastrowid

def add_cutoff(mdb,event):
	sql = '''   INSERT OR IGNORE INTO cutoffs (obsidtxt,obsid,cz,binning, method, cutb0Q0, cutb0Q1, cutb0Q2, cutb0Q3,cutb1Q0, cutb1Q1, cutb1Q2, cutb1Q3, cutb2Q0, cutb2Q1, cutb2Q2, cutb2Q3)VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'''
	cur=mdb.cursor()
	cur.execute(sql,event)
	return cur.lastrowid

def add_cutoff_veto(mdb,event):
	sql = '''   INSERT OR IGNORE INTO cutoffs (obsidtxt, obsid,cz,binning, method, cutb0Q0, cutb0Q1, cutb0Q2, cutb0Q3) VALUES(?,?,?,?,?,?,?,?,?)'''
	cur=mdb.cursor()
	cur.execute(sql,event)
	return cur.lastrowid

def saa_times(saa_start, saa_end, temp_startbin, temp_stopbin, bins_time, peakind):
	if len(saa_start) >= 1:
		# if(len(temp_startbin)==0) and (len(temp_stopbin)==0):
		# 	pre_edge=bins_time[0]
		# 	post_edge=bins_time[-1]
		# elif len(temp_startbin)==0:
		# 	pre_edge=bins_time[0]
		# elif len(temp_stopbin)==0:
		# 	post_edge=bins_time[-1]
		# else:
			# print("Peak and SAA times:" , bins_time[peakind], saa_start, saa_end)
			
			for i in range(len(saa_start)):
				if bins_time[peakind] <= saa_start[0]:
					pre_edge = np.abs(bins_time[peakind] - saa_start[0])
					post_edge = -1
					print(bins_time[peakind], saa_start[0])
					break
				elif bins_time[peakind] >= saa_start[i] and bins_time[peakind] <= saa_end[i]:
					pre_edge = 0
					post_edge = 0
					print(saa_start[i], bins_time[peakind], saa_end[i])
					break
				elif i+1< len(saa_start) and (bins_time[peakind] >= saa_end[i] and bins_time[peakind] <= saa_start[i+1]):
					pre_edge = np.abs(bins_time[peakind] - saa_start[i+1])
					post_edge = np.abs(bins_time[peakind] - saa_end[i])
					print(saa_end[i], bins_time[peakind], saa_start[i+1])
					break
				else:
					pre_edge = -1
					post_edge = np.abs(bins_time[peakind] - saa_end[i])
					print(saa_end[i], bins_time[peakind])
	else:
		pre_edge = -1
		post_edge = -1
	
	return post_edge, pre_edge

def grb_search(args, outpaths, orbitinfo, startbins, stopbins, dirname, method, saa_start, saa_end, thresh):
	# ## making table for saving detected peaks
	orbit = outpaths[0].split('/')[-1]
	numgrbs = 0
	for c_tbin, n_tbin in enumerate(args.tbin):
		peak_maps = []
		bins_time = []
		time,lc_all,mask_all = [],[],[]
		cutoffs_tosave = []
		counts = []
		for band in range(3):
			tbl = Table.read('/home/ashwin/False_Alarm_Rate1/data/evt_files/'+str(orbit)+'/'+str(orbit)+'_'+str(n_tbin)+'_'+str(band)+'_Q'+str(0)+'_detrended.fits')
			y = np.shape(np.array(tbl['countrate']))
			counts.append(y)
		min_yc = np.min(counts)
		for band in range(3):
			cutoff_bin = [0,0,0,0]
			lc_band,mask_band = [],[] ##loading quadrant lcs and masks for each band
			for quad in range(4):
				tbl = Table.read('/home/ashwin/False_Alarm_Rate1/data/evt_files/'+str(orbit)+'/'+str(orbit)+'_'+str(n_tbin)+'_'+str(band)+'_Q'+str(quad)+'_detrended.fits')
				lc_band.append(np.array(tbl['countrate'])[:min_yc])
				mask_band.append(np.array(tbl['mask'])[:min_yc])
				time = np.array(tbl['time'])[:min_yc]
			lc_all.append(lc_band)
			mask_all.append(mask_band)
			bins_time = (time[1:]+time[:-1])*0.5
			if(method == 'cumsum'):
				##histogram for cutoff calculation
				hist_bin = preproc.gethist(lc_band,mask_band)
				cutoff_bin = preproc.getcutoffrate1(n_tbin, hist_bin, args.far,args.timespan,lc_band,mask_band)
				cutoffs_tosave.append(cutoff_bin)
				##get peak map
				peak_ind_bin, peak_map = preproc.getpeaks(lc_band, cutoff_bin)
				peak_maps.append(peak_map)
			elif(method == 'nsigma'):
				cutoff_bin = preproc.getcutoffrate2(lc_band,mask_band, thresh)
				peak_ind_bin, peak_map = preproc.getpeaks(lc_band, cutoff_bin)
				peak_maps.append(peak_map)
				cutoffs_tosave.append(cutoff_bin)
				# print(f"here: {np.sum(peak_map, axis=1)} crossings")
				
		# cutoffs_tosave = cutoffs_tosave.astype(float)
		obsidtxt=dirname+'_'+str(n_tbin)+'_'+method+'_czti'
		evs=[obsidtxt, dirname, 'czti',n_tbin, method+'_czti']+cutoffs_tosave[0]+cutoffs_tosave[1]+cutoffs_tosave[2]

		# lastrid=add_cutoff(mdb,evs)
		# lastrid=add_cutoff(bsl,evs)
		lc_all = np.array(lc_all)
		# print("Doing CZTI analysis, binnning time : "+str(args.tbin[c_tbin])+", method : "+method)
		band_peak_maps=np.vstack([peak_maps[0],peak_maps[1],peak_maps[2]])
		# print("Shape of band peak map: ",band_peak_maps.shape)
		peak_ind_bin = np.where(np.sum(band_peak_maps, axis=0) >= 4.0)
		# print("peak_ind_bin:",peak_ind_bin)
		for orb in range(len(outpaths)):
			global start_orbit_czti, end_orbit_czti
			start_orbit_czti = orbitinfo['start_of_orbit'][orb]
			end_orbit_czti = orbitinfo['end_of_orbit'][orb]
			startbins = np.array(startbins)
			stopbins = np.array(stopbins)
			# print("startbins, stopbins", startbins,stopbins)
		storage = []
		if peak_ind_bin[0].tolist():
			peakinds = peak_ind_bin[0]
			# print('peakinds',peakinds)
			total_counts = preproc.total_countrate(bins_time,peakinds,peak_ind_bin[0],lc_all)
			# print('total counts:', total_counts)
			print("GRB_time   Rank   Detected_in_quads  rate")
			print('------------------------------------------')
			det_quad = []
			for counter, peakind in enumerate(peakinds):
				if(peakind == len(bins_time)):
					continue
				quads = np.where(band_peak_maps[:, peakind])
				det_quad.append(quads)
				grb_orbit=''
				for orb in range(len(outpaths)):
					start_orbit = orbitinfo['start_of_orbit'][orb]
					end_orbit = orbitinfo['end_of_orbit'][orb]
					if(bins_time[peakind]>=start_orbit and bins_time[peakind]<=end_orbit):
						grb_orbit = outpaths[orb]
				orbitname = os.path.basename(grb_orbit)
				tplot,lc,mask = preproc.grb_lc(lc_all, mask_all, time, peakind, n_tbin)
				stat_significance,band_significance = preproc.get_significance(lc,mask,total_counts[counter])
				stat_significance = np.round(stat_significance,1)
				band_significance = np.round(band_significance,1)
				#Find detected quadrants and band and get rank
				rank,numquads0,numquads1,numquads2 = preproc.get_rank(quads[0])
				numquads = np.unique((numquads0+numquads1+numquads2))
				print((str(round(bins_time[peakind],0))+'   '+str(rank)+'   '+','.join(str(f) for f in numquads)+'  '+str(np.sum(lc_all[:,:,peakind]))))
				#Find distance from SAA
				startbins = np.array(startbins)
				stopbins = np.array(stopbins)
				temp_startbin = startbins[startbins<=bins_time[peakind]]
				temp_stopbin = stopbins[stopbins>=bins_time[peakind]]
				lc_all = np.round(lc_all,2)

				post_edge, pre_edge = saa_times(saa_start, saa_end, temp_startbin, temp_stopbin, bins_time, peakind)

				print("Distance from SAA", post_edge)
				print("Distance to SAA", pre_edge)

				
				method_tag = {'cumsum' : 'C', 'nsigma' : 'N'}
				eventid = 'EZ'+method_tag[method]+'{0:1.1e}'.format(n_tbin)+'_'+str(round(bins_time[peakind],0))
				event = [eventid, round(bins_time[peakind],0), n_tbin, method, orbitname, dirname, round(post_edge,0), round(pre_edge,0), rank,
						','.join(str(f) for f in numquads0), ','.join(str(f) for f in numquads1),','.join(str(f) for f in numquads2),
						band_significance[0], band_significance[1], band_significance[2],
						stat_significance[0],stat_significance[1],stat_significance[2],stat_significance[3],
						stat_significance[4],stat_significance[5],stat_significance[6],stat_significance[7],
						stat_significance[8],stat_significance[9],stat_significance[10],stat_significance[11],
						lc_all[0][0][peakind],lc_all[0][1][peakind],lc_all[0][2][peakind],lc_all[0][3][peakind],
						lc_all[1][0][peakind],lc_all[1][1][peakind],lc_all[1][2][peakind],lc_all[1][3][peakind],
						lc_all[2][0][peakind],lc_all[2][1][peakind],lc_all[2][2][peakind],lc_all[2][3][peakind]]
				storage.append(event)
				putInBS = False
				if method=='nsigma' and n_tbin == 0.1:
					# Need to check for readout dips in nsigma and 0.1 tbins. 
					print(glob.glob('/home/ashwin/False_Alarm_Rate1/data/evt_files/'+orbit+'/*'+orbit+'*quad_clean.evt'))
					eventfile = glob.glob('/home/ashwin/False_Alarm_Rate1/data/evt_files/'+orbit+'/*'+orbit+'*quad_clean.evt')[0]
					hdu = fits.open(eventfile)
					readouts = [False, False, False, False] # Whether we have readout dip in each quadrant
					for quad in range(4):
						data = hdu[quad+1].data
						t_min = min(data['time'])  
						t_max = max(data['time']) 
						emin=20.
						emax=200.
						ebin=10.
						ebins = np.arange(emin, emax+ebin/100., ebin)
						times = np.arange(t_min, t_max+n_tbin, n_tbin)
						h = np.histogram2d(data['time'], data['energy'], bins=(times, ebins))[0]
						h = np.nansum(h, axis=1) / n_tbin # At this point h is the light curve
						ind = int((bins_time[peakind]-t_min)/n_tbin) # Find which index to use, this calculation is required since bins_time[0] and t_min don't match.
						for i in range(-5, 5):
							# Check in +- 5 bins
							if h[ind-i] == 0 and h[ind-i+1] == 0:
								# Readout is when we have 2-3 consecutive 0s. We only check for 2 consecutive 0s
								readouts[quad] = True
								break
						# Check for readout dip in > 2 quadrants
					if sum(readouts) > 3:
						# Put this in bs
						putInBS = True
				if pre_edge == 0:
					# lastrowid = add_czti_event(bsl, event)
					print("In SAA")
				elif putInBS:
					# lastrowid = add_czti_event(bsl, event)
					print("Readout")
				else:
					# lastrowid = add_czti_event(mdb, event)
					print('Good Event')
					numgrbs = numgrbs+1
		else:
			numgrbs = numgrbs + 0
	return numgrbs, storage

def grb_veto_search(mdb, bsl, args, outpaths, orbitinfo, startbins, stopbins, dirname, method, saa_start, saa_end):
	numgrbs_veto = 0
	for c_tbin,n_tbin in enumerate(args.tbin):
		if(n_tbin == 0.1):
			continue
		bins_time = []
		cutoff_bin = [0,0,0,0]
		time,lc_all,mask_all = [],[],[]
		try:
			for quad in range(4):
				tbl = Table.read('data/local_level2/'+dirname+'/czti/combined_veto_'+str(n_tbin)+'_Q'+str(quad)+'_detrended.fits')
				lc_all.append(np.array(tbl['vetocounts']))
				mask_all.append(np.array(tbl['vetomask']))
				time = np.array(tbl['time'])
		except:
			print("Error in directory "+dirname+" and tbin "+str(n_tbin)+" in veto")
			continue
			##histogram for cutoff calculation
		bins_time = (time[1:]+time[:-1])*0.5
		if(method == 'cumsum'):
			##histogram for cutoff calculation
			hist_bin = preproc.gethist(lc_all,mask_all)
			cutoff_bin = preproc.getcutoffrate1(n_tbin, hist_bin, args.far,args.timespan,lc_all,mask_all)
			# for quad in range(4):
			# 	bins = hist_bin[quad][0]
			# 	bins = 0.5*(bins[1:]+bins[:-1])
			# 	fig1 = plt.figure()
			# 	plt.plot(bins,hist_bin[quad][1],label='Quad '+str(quad))
			# 	plt.axvline(cutoff_bin[quad],color='black',linestyle='dashed',label='Cutoff rate')
			# 	plt.title(str(band))
			# 	plt.xlabel('Rates')
			# 	plt.ylabel('Rate histogram')
			# 	plt.savefig(dirname+'_'+str(n_tbin)+'_'+str(band)+'_'+str(quad)+'_histograms.png')
			# 	plt.close()
			##get peak map
			peak_ind_bin, peak_map = preproc.getpeaks(lc_all, cutoff_bin)
		elif(method == 'nsigma'):
			cutoff_bin = preproc.getcutoffrate2(lc_all,mask_all, 5)
			peak_ind_bin, peak_map = preproc.getpeaks(lc_all, cutoff_bin)
		obsidtxt=dirname+'_'+str(n_tbin)+'_'+method+'_veto'
		evs=[obsidtxt, dirname, 'veto',n_tbin, method+'_veto', cutoff_bin[0],cutoff_bin[1],cutoff_bin[2],cutoff_bin[3] ]
		lastrid=add_cutoff_veto(mdb,evs)
		lastrid=add_cutoff_veto(bsl,evs)
		lc_all = np.array(lc_all)
		lc = np.where(lc_all <= 0)
		print("Doing VETO analysis, binnning time : "+str(args.tbin[c_tbin])+", method : "+method)
		peak_ind_bin = np.where(np.sum(peak_map, axis=0) >= 2.0)
		if peak_ind_bin[0].tolist():			
			peakinds = peak_ind_bin[0]
			for counter, peakind in enumerate(peakinds):
				start_orbit = start_orbit_czti
				end_orbit = end_orbit_czti
				if(bins_time[peakind]>=start_orbit and bins_time[peakind]<=end_orbit):
					total_counts = preproc.total_countrate_veto(bins_time,peakinds,peak_ind_bin[0],lc_all)
					print("GRB_time   Significance   Detected_in_quads   rate")
					print('---------------------------------------------------')
					det_quad = []
					if(peakind == len(bins_time)):
						continue
					quads = np.where(peak_map[:, peakind])
					det_quad.append(quads)
					grb_orbit = ''
					for orb in range(len(outpaths)):
						# start_orbit = orbitinfo['veto_start'][orb]
						# end_orbit = orbitinfo['veto_end'][orb]
						if(bins_time[peakind]>=start_orbit and bins_time[peakind]<=end_orbit):
							grb_orbit = outpaths[orb]
					orbitname = os.path.basename(grb_orbit)
					tplot,lc,mask = preproc.grb_lc_veto(lc_all,mask_all,time,peakind,n_tbin)
					stat_significance,tot_significance = preproc.get_significance_veto(lc,mask,total_counts[counter])
					print((str(round(bins_time[peakind],0))+'   '+str(tot_significance)+'   '+','.join(str(f) for f in quads[0])+'  '+str(np.sum(lc_all[:,peakind]))))
					stat_significance = np.round(stat_significance,1)
					tot_significance = np.round(tot_significance,1)

					#Find distance from SAA
					starting_bins = np.array(startbins)
					stopping_bins = np.array(stopbins)
					startbins = np.array(start_orbit_czti)
					stopbins = np.array(end_orbit_czti)
					temp_startbin = startbins[startbins<=bins_time[peakind]]
					temp_stopbin = stopbins[stopbins>=bins_time[peakind]]
					lc_all = np.round(lc_all,2)

					post_edge, pre_edge = saa_times(saa_start, saa_end, temp_startbin, temp_stopbin, bins_time, peakind)
					print("Distance from SAA", post_edge)
					print("Distance to SAA", pre_edge)
					##################plot grb function here
					method_tag = {'cumsum' : 'C', 'nsigma' : 'N'}
					eventid = 'EV'+method_tag[method]+'{0:1.1e}'.format(n_tbin)+'_'+str(round(bins_time[peakind],0))
					event = [eventid, round(bins_time[peakind],0), n_tbin, method+'_veto', orbitname, dirname, round(post_edge,0), round(pre_edge,0),
							','.join(str(f) for f in quads[0]), tot_significance, 
							stat_significance[0], stat_significance[1],stat_significance[2],stat_significance[3], 
							lc_all[0][peakind],lc_all[1][peakind],lc_all[2][peakind],lc_all[3][peakind]]
					if pre_edge == 0:
						lastrowid = add_veto_event(bsl, event)
						print("In SAA, added to black_sheep")
					else:
						lastrowid = add_veto_event(mdb, event)
						print('Good Event, added to motherload')
						numgrbs_veto = numgrbs_veto+1
				else:
					numgrbs_veto = numgrbs_veto + 0
	return numgrbs_veto	

def topN_search(mdb, bsl, args, outpaths, orbitinfo, startbins, stopbins, dirname, saa_start, saa_end):
	## making table for saving detected peaks
	numgrbs = 0
	for o, orbit in enumerate(outpaths):
		orbitname = os.path.basename(orbit)
		for c_tbin,n_tbin in enumerate(args.tbin):
			peak_maps = []
			bins_time = []
			cutoff_bin = [0,0,0,0]
			time,lc_all,mask_all = [],[],[]
			try:
				for band in range(3):
					lc_band,mask_band = [],[] ##loading quadrant lcs and masks for each band
					for quad in range(4):
						tbl = Table.read(orbit+'/'+orbitname+'_'+str(n_tbin)+'_'+str(band)+'_Q'+str(quad)+'_detrended.fits')
						lc_band.append(np.array(tbl['countrate']))
						mask_band.append(np.array(tbl['mask']))
						time = np.array(tbl['time'])
					lc_all.append(lc_band)
					mask_all.append(mask_band)
					bins_time = (time[1:]+time[:-1])*0.5
					##get peak map
					peak_ind_bin, peak_map = preproc.getpeakstopN(lc_band,3)
					peak_maps.append(peak_map)
				lc_all = np.array(lc_all)
				print("Doing CZTI analysis, binnning time : "+str(args.tbin[c_tbin])+", method : topN, orbit : "+orbitname)
				band_peak_maps=np.vstack([peak_maps[0],peak_maps[1],peak_maps[2]])
				peak_ind_bin = np.where(np.sum(band_peak_maps, axis=0) >= 4.0)

				if peak_ind_bin[0].tolist():
					peakinds = peak_ind_bin[0]
					total_counts = preproc.total_countrate(bins_time,peakinds,peak_ind_bin[0],lc_all)
					print("GRB_time   Rank   Detected_in_quads  rate")
					print('------------------------------------------')
					det_quad = []
					for counter, peakind in enumerate(peakinds):
						if(peakind == len(bins_time)):
							continue
						quads = np.where(band_peak_maps[:, peakind])
						det_quad.append(quads)
						grb_orbit=orbit
						orbitname = os.path.basename(grb_orbit)
						tplot,lc,mask = preproc.grb_lc(lc_all, mask_all,time,peakind,n_tbin)
						stat_significance,band_significance = preproc.get_significance(lc,mask,total_counts[counter])
						stat_significance = np.round(stat_significance,1)
						band_significance = np.round(band_significance,1)
						#Find detected quadrants and band and get rank
						rank,numquads0,numquads1,numquads2 = preproc.get_rank(quads[0])
						numquads = np.unique((numquads0+numquads1+numquads2))
						print((str(round(bins_time[peakind],0))+'   '+str(rank)+'   '+','.join(str(f) for f in numquads)+'  '+str(np.sum(lc_all[:,:,peakind]))))
						#Find distance from SAA
						startbins = np.array(startbins)
						stopbins = np.array(stopbins)
						temp_startbin = startbins[startbins<=bins_time[peakind]]
						temp_stopbin = stopbins[stopbins>=bins_time[peakind]]
						lc_all = np.round(lc_all,2)

						post_edge, pre_edge = saa_times(saa_start, saa_end, temp_startbin, temp_stopbin, bins_time, peakind)

						print("Distance from SAA", post_edge)
						print("Distance to SAA", pre_edge)
						eventid = 'EZT'+'{0:1.1e}'.format(n_tbin)+'_'+str(round(bins_time[peakind],0))
						event = [eventid, round(bins_time[peakind],0), n_tbin, 'topN', orbitname, dirname, round(post_edge,0), round(pre_edge,0), rank,
								','.join(str(f) for f in numquads0), ','.join(str(f) for f in numquads1),','.join(str(f) for f in numquads2),
								band_significance[0], band_significance[1], band_significance[2],
								stat_significance[0],stat_significance[1],stat_significance[2],stat_significance[3],
								stat_significance[4],stat_significance[5],stat_significance[6],stat_significance[7],
								stat_significance[8],stat_significance[9],stat_significance[10],stat_significance[11],
								lc_all[0][0][peakind],lc_all[0][1][peakind],lc_all[0][2][peakind],lc_all[0][3][peakind],
								lc_all[1][0][peakind],lc_all[1][1][peakind],lc_all[1][2][peakind],lc_all[1][3][peakind],
								lc_all[2][0][peakind],lc_all[2][1][peakind],lc_all[2][2][peakind],lc_all[2][3][peakind]]

						if pre_edge == 0:
							lastrowid = add_czti_event(bsl, event)
							print("In SAA")
						else:
							lastrowid = add_czti_event(mdb, event)
							print('Good Event')
							numgrbs = numgrbs+1
						
				else:
					numgrbs = numgrbs + 0
			except:
				continue
	return numgrbs


def topN_veto_search(mdb, bsl, args, outpaths, orbitinfo, startbins, stopbins, dirname, saa_start, saa_end):
	## making table for saving detected peaks
	numgrbs_veto = 0
	for o,orbit in enumerate(outpaths):
		orbitname = os.path.basename(orbit)
		for c_tbin,n_tbin in enumerate(args.tbin):
			if(n_tbin == 0.1):
				continue
			bins_time = []
			cutoff_bin = [0,0,0,0]
			time,lc_all,mask_all = [],[],[]
			try:
				for quad in range(4):
					tbl = Table.read(orbit+'/'+orbitname+'_veto_'+str(n_tbin)+'_Q'+str(quad)+'_detrended.fits')
					lc_all.append(np.array(tbl['vetocounts']))
					mask_all.append(np.array(tbl['vetomask']))
					time = np.array(tbl['time'])
			except:
				print(("Error : No detrended fits file found for veto orbit "+orbitname+" of directory "+
					dirname+" tbin and quad "+str(n_tbin)+" , "+str(quad)))
				continue
			bins_time = (time[1:]+time[:-1])*0.5
			##get peak map
			peak_ind_bin, peak_map = preproc.getpeakstopN(lc_all,3)
			lc_all = np.array(lc_all)
			print("Doing VETO analysis, binnning time : "+str(args.tbin[c_tbin])+", method : topN, orbit : "+orbitname)
			peak_ind_bin = np.where(np.sum(peak_map, axis=0) >= 2.0)
			# np.save('arr_'+dirname+'_'+orbitname,peak_map)
			# np.save('arr_time_'+dirname+'_'+orbitname,time)
			if peak_ind_bin[0].tolist():
				peakinds = peak_ind_bin[0]
				total_counts = preproc.total_countrate_veto(bins_time,peakinds,peak_ind_bin[0],lc_all)
				for counter, peakind in enumerate(peakinds):
					start_orbit = start_orbit_czti
					end_orbit = end_orbit_czti
					if(bins_time[peakind]>=start_orbit and bins_time[peakind]<=end_orbit):
						print("GRB_time   Significance   Detected_in_quads   rate")
						print('---------------------------------------------------')
						det_quad = []
						if(peakind == len(bins_time)):
							continue
						quads = np.where(peak_map[:, peakind])
						det_quad.append(quads)
						grb_orbit=orbit
						orbitname = os.path.basename(grb_orbit)
						tplot,lc,mask = preproc.grb_lc_veto(lc_all,mask_all,time,peakind,n_tbin)
						stat_significance,tot_significance = preproc.get_significance_veto(lc,mask,total_counts[counter])
						print((str(round(bins_time[peakind],0))+'   '+str(tot_significance)+'   '+','.join(str(f) for f in quads[0])+'  '+str(np.sum(lc_all[:,peakind]))))
						stat_significance = np.round(stat_significance,1)
						tot_significance = np.round(tot_significance,1)

						#Find distance from SAA
						startbins = np.array(startbins)
						stopbins = np.array(stopbins)
						temp_startbin = startbins[startbins<=bins_time[peakind]]
						temp_stopbin = stopbins[stopbins>=bins_time[peakind]]
						lc_all = np.round(lc_all,2)

						post_edge, pre_edge = saa_times(saa_start, saa_end, temp_startbin, temp_stopbin, bins_time, peakind)

						print("Distance from SAA", post_edge)
						print("Distance to SAA", pre_edge)
						eventid = 'EVT'+'{0:1.1e}'.format(n_tbin)+'_'+str(round(bins_time[peakind],0))
						event = [eventid, round(bins_time[peakind],0), n_tbin, 'topN_veto', orbitname, dirname, round(post_edge,0), round(pre_edge,0),
								','.join(str(f) for f in quads[0]), tot_significance, 
								stat_significance[0], stat_significance[1],stat_significance[2],stat_significance[3], 
								lc_all[0][peakind],lc_all[1][peakind],lc_all[2][peakind],lc_all[3][peakind]]
						if pre_edge == 0:
							lastrowid = add_veto_event(bsl, event)
							print("In SAA")
						else:
							lastrowid = add_veto_event(mdb, event)
							print('Good Event')
							numgrbs_veto = numgrbs_veto+1

					else:
						numgrbs_veto = numgrbs_veto + 0
	return numgrbs_veto
