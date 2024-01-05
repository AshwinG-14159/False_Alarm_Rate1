import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import glob
import sqlite3
import t90_mod as calct90
from matplotlib.ticker import MultipleLocator,AutoMinorLocator
from PyPDF2 import PdfReader, PdfWriter

def calc_stats(lc, times, tmark, tbin, t_min, t_max, plotfilename, eventfile, plotfile):
	times = np.array(times, dtype=object)
	lc = np.array(lc, dtype=object)
	try:
		lcnew = np.sum(lc,axis=0)
		timesnew = times[0]
	except:
		minlen = np.min([len(lc[0]),len(lc[1]),len(lc[2]),len(lc[3])])
		lc = np.array([lc[0][0:minlen],lc[1][0:minlen],lc[2][0:minlen],lc[3][0:minlen]])
		lcnew = np.sum(lc,axis=0)
		timesnew = times[0][0:minlen]
	tran_start = tmark[0]-tbin
	tran_end = tmark[-1]+tbin
	quads = [0,1,2,3]
	# print(len(lcnew), len(timesnew))
	#print plotfilename
	#print len(lc), len(times), tmark, t_min, t_max
	t90, bg, rate, counts = calct90.iterate_t90(eventfile, plotfilename, t_min, t_max, tran_start, tran_end, quads, tbin, lcnew, timesnew, str(tmark[0]), plotfile)
	return t90, bg, rate, counts

def make_plot(plotfile,h,e_bin,tbin,quad,elower,eupper,ebins,tmin,tmax,times,tmark, unit="Counts", filename="lc", title="Lightcurve", specscale="linear"):
	lc = np.nansum(h, axis=1) / tbin
	spec = np.nansum(h, axis=0) / e_bin
	fig1 = plt.figure()
	gs = gridspec.GridSpec(2, 2, width_ratios=[3,1], height_ratios=[2,1],hspace=0.2,wspace=0.05)
	
	# Energy-time plot
	ax_et = fig1.add_subplot(gs[0])
	if('veto' in filename):
		ax_et.imshow(np.transpose(h), extent=(float(tmin), float(tmax), elower, eupper), origin='lower', aspect='auto',
					vmin = np.min(h[h>=0]), vmax = np.max(h[:,0:100]))  # left, right, bottom, top
	else:
		ax_et.imshow(np.transpose(h), extent=(float(tmin), float(tmax), elower, eupper), origin='lower', aspect='auto')  # left, right, bottom, top		
	fig1.suptitle("Quadrant " + ['A', 'B', 'C', 'D'][quad])
	ax_et.set_ylabel('Energy (keV)', color='black',fontsize='small')
	#ax_et.set_xlabel('GPS Time (sec)', color='black')
	ax_et.xaxis.set_label_position('top')
	for tm in tmark:
		ax_et.annotate('',fontsize='xx-small',xy=(float(tm),eupper),xytext=(float(tm),eupper+15),arrowprops=dict(facecolor='black',shrink=0.001))
		ax_et.annotate('',fontsize='xx-small',xy=(float(tm),elower),xytext=(float(tm),elower-15),arrowprops=dict(facecolor='black',shrink=0.001))
	#	ax_et.axvline(float(tm), color='black', linestyle='dashed', label='Marked times')
	ax_et.set_ylim(elower, eupper)
	ax_et.yaxis.set_minor_locator(AutoMinorLocator(5))
	ax_et.xaxis.set_minor_locator(MultipleLocator(tbin*2.5))
	offset=ax_et.xaxis.get_offset_text()
	offset.set_size('x-small')
	#ax_et.ticklabel_format(style='plain')
	plt.xticks(fontsize='x-small')
	plt.yticks(fontsize='x-small')		
	
	# Energy spectrum
	ax_e = fig1.add_subplot(gs[1])
	eplot = 0.5 * (ebins[1:] + ebins[:-1])
	#eplot = np.arange(0,256,1)
	#print len(np.nansum(h,axis=0)/e_bin)
	ax_e.plot(spec, eplot,'o-',fillstyle='none',markersize=3,lw=1,color='black')
	ax_e.set_xscale(specscale)
	ax_e.set_ylim(elower, eupper)
	ax_e.xaxis.set_ticks_position('top') 
	plt.setp(ax_e.xaxis.get_majorticklabels(), rotation=45 )
	ax_e.set_xlabel(unit+"/keV", color='black',fontsize='small')
	#ax_e.set_ylabel('Energy (keV)', color='black',fontsize='small')
	ax_e.yaxis.set_ticks_position('right')
	ax_e.yaxis.set_minor_locator(AutoMinorLocator(5))
	#offset=ax_e.xaxis.get_offset_text()
	# #offset.set_size('x-small')
	plt.xticks(fontsize='x-small')
	plt.yticks(fontsize='x-small')
	ax_e.grid(alpha=0.3)	
	
	# Time (lightcurve)
	ax_t = plt.subplot(gs[2])
	tplot = times
	if(len(tplot)!=len(lc)):
		tplot = np.array(0.5 * (times[1:] + times[:-1]))
	ax_t.plot(tplot, lc,lw=1,color='black')
	ax_t.set_ylabel(unit+"/sec", color='black',fontsize='small')
	ax_t.set_xlabel("Trigger Time (sec)", color='black',fontsize='small')
	ylims = plt.ylim()
	yupper = ylims[1] + 0.30 * (ylims[1] - ylims[0])
	for tm in tmark:
	        ax_t.axvline(float(tm), color='black', alpha=0.8,lw=0.8, linestyle='dashed', label='Marked times')
	ax_t.set_ylim(ylims[0], yupper)
	ax_t.set_xlim(tmin, tmax)
	ax_t.legend(fontsize='x-small', framealpha=0.5)
	ylims = ax_t.get_ylim()
	ax_t.yaxis.set_minor_locator(AutoMinorLocator(5))
	ax_t.xaxis.set_minor_locator(MultipleLocator(tbin*2.5))
	offset=ax_t.xaxis.get_offset_text()
	offset.set_size('x-small')
	plt.xticks(fontsize='x-small')
	plt.yticks(fontsize='x-small')
	ax_t.grid(alpha=0.3)
		
	# Count rate histogram
	ax_stats = fig1.add_subplot(gs[3])
	lc_stats = np.histogram(lc, bins=np.linspace(ylims[0], ylims[1], 20))
	plot_ticks = 0.5 * (lc_stats[1][1:] + lc_stats[1][:-1])
	ax_stats.plot(lc_stats[0], plot_ticks,lw=1,color='black')
	ax_stats.set_xlabel("number of bins", color='black',fontsize='small')
	#ax_stats.set_ylabel(unit+"/sec", color='black',fontsize='small')
	#ax_stats.yaxis.set_label_position('right')
	ax_stats.yaxis.set_ticks_position('right')
	ax_stats.yaxis.set_minor_locator(AutoMinorLocator(5))
	ax_stats.xaxis.set_minor_locator(AutoMinorLocator(5))
	plt.xticks(fontsize='x-small')
	plt.yticks(fontsize='x-small')
	ax_stats.grid(alpha=0.3)
		
	plotfile.savefig()
	plt.close()
	return lc

def plot_czti_lc(plotflag, statfile, pdffolder, dirname, orbitname, tmark, tbin, tmin, tmax, emin=20., emax=200., ebin=10.):
	tmark = tmark.astype(float)
	tmin = float(tmin)
	tmax = float(tmax)
	tbin = float(tbin)
	plotfilename = "{stem}_lightcurves_czti.pdf".format(stem=pdffolder+'/'+str(tmark[0])+'_'+str(tbin))
	if(len(glob.glob(plotfilename)) != 0 and plotflag==False):
		return Table.read(pdffolder+'/statfile.csv')
	plotfile = PdfPages(plotfilename)
	eventfile = glob.glob('/home/czti/grb_search/data/local_level2/'+dirname+'/czti/orbit/'+orbitname+'/*quad_clean.evt')[0]
	hdu = fits.open(eventfile)
	tmin = tmin - 50.0*tbin
	tmax = tmax + 50.0*tbin
	ebins = np.arange(emin, emax+ebin/100., ebin)
	num_ebins = len(ebins)-1
	lc_all = []
	times_all = []
	for quad in range(4):
		quadname = ['A', 'B', 'C', 'D'][quad]
		data = hdu[quad+1].data
		t_min = max((min(data['time']), tmin) )  # Start at datastart or args.tmin, whichever comes LATER
		t_max = min((max(data['time']), tmax) )  # End at dataend or args.tmax, whichever comes FIRST
		times = np.arange(t_min, t_max+tbin, tbin)
		times_all.append(times)
		if(len(times) == 0):
			statf_row = [orbitname,'czti',tbin,0,0,0,-100]
			statfile.add_row(statf_row)
			plotfile.close()
			return statfile
		h = np.histogram2d(data['time'], data['energy'], bins=(times, ebins))
		h = h[0]
		# print(t_min, t_max, tbin)
		# print("czti", len(h),len(times))
		sigmas = np.zeros(num_ebins)
		medians = np.zeros(num_ebins)
		hnorm = np.zeros(np.shape(h))
		hsub = np.zeros(np.shape(h))
		for count in range(num_ebins):
			calc = sigma_clipped_stats(h[:,count], sigma=3.0, maxiters=3)
			medians[count] = calc[1]
			if calc[2] != 0:
				sigmas[count] = calc[2]
			else:
				sigmas[count] = 1.0
			hsub[:,count] = h[:,count] - medians[count]
			hnorm[:,count] = (h[:,count] - medians[count]) / sigmas[count]
		lc = make_plot(plotfile,h,ebin,tbin,quad,emin,emax,ebins,t_min,t_max,times,tmark, filename="{q}_lc_czti".format(q=quadname), title="LC:Quad {q}".format(q=quadname), specscale='log')
		lcsub = make_plot(plotfile,hsub,ebin,tbin,quad,emin,emax,ebins,t_min,t_max,times,tmark, unit="Excess counts", filename="{q}_lcsub_czti".format(q=quadname), title="LC-median subtracted:Quad {q}".format(q=quadname))
		lcnorm = make_plot(plotfile,hnorm,ebin,tbin,quad,emin,emax,ebins,t_min,t_max,times,tmark, unit="Strength", filename="{q}_lcnorm_czti".format(q=quadname), title="LC-normalised:{q}".format(q=quadname))
		lc_all.append(lc)
		# Combined lightcurves
		fig = plt.figure()
		gs = gridspec.GridSpec(2, 1, height_ratios=[2,1],hspace=0.2)
		ax_top = fig.add_subplot(gs[0])
		ax_top.set_title("Lightcurve:"+" Quadrant " + quadname)
		tplot = 0.5 * (times[1:] + times[:-1])
		ax_top.plot(tplot, lc, alpha=0.9, label="{emin:0.1f}-{emax:0.1f} keV LC".format(emin=emin, emax=emax))
		ax_top.plot(tplot, lcsub, alpha=0.9, label="Median spectrum subtracted LC")
		ax_top.set_ylabel("Counts/sec", color='black')
		#ax_top.set_xlabel("Time (sec)", color='red')
		ylims = ax_top.get_ylim()
		yupper = ylims[1] + 0.20 * (ylims[1] - ylims[0])
		for tm in tmark:
			ax_top.axvline(tm, color='black', alpha=0.8, lw=0.8, linestyle='dashed', label='Marked times')
		ax_top.set_ylim(ylims[0], yupper)
		ax_top.set_xlim(t_min, t_max)
		plt.xticks(fontsize='small')
		plt.yticks(fontsize='small')
		ax_top.xaxis.set_minor_locator(MultipleLocator(tbin*2.5))	
		ax_top.yaxis.set_minor_locator(AutoMinorLocator(5))
		offset=ax_top.xaxis.get_offset_text()
		offset.set_size('x-small')
		ax_top.legend(fontsize='x-small', handlelength=1, framealpha=0.5)
		ax_top.grid(alpha=0.3)
		#ax_bottom = plt.subplot(gs[1])
		ax_bottom=fig.add_subplot(gs[1],sharex=ax_top)
		ax_bottom.plot(tplot, lcnorm, label="Mean spectrum subtracted + normalized LC")
		ax_bottom.set_xlabel("Trigger Time (sec)", color='black')
		ax_bottom.set_ylabel("Significance")
		ylims = ax_bottom.get_ylim()
		yupper = ylims[1] + 0.40 * (ylims[1] - ylims[0])
		for tm in tmark:
			ax_bottom.axvline(tm, color='black', alpha=0.8,lw=0.8,linestyle='dashed', label='Marked times')
		ax_bottom.set_ylim(ylims[0], yupper)
		ax_bottom.set_xlim(t_min, t_max)
		plt.xticks(fontsize='small')
		plt.yticks(fontsize='small')
		ax_bottom.xaxis.set_minor_locator(MultipleLocator(tbin*2.5))	
		ax_bottom.yaxis.set_minor_locator(AutoMinorLocator(5))
		offset=ax_bottom.xaxis.get_offset_text()
		offset.set_size('x-small')
		ax_bottom.legend(fontsize='x-small', handlelength=1, framealpha=0.5)
		ax_bottom.grid(alpha=0.3)
		plotfile.savefig()
		plt.close()
	t90, bg, rate, counts = calc_stats(lc_all, times_all, tmark, tbin, t_min, t_max, plotfilename, eventfile, plotfile)
	statf_row = [orbitname,'czti',tbin,bg,rate,counts,t90]
	statfile.add_row(statf_row)
	fig = plt.figure()
	for quad in range(4):
		# plot all 4 lightcurves on a single plot
		quadname = ['A', 'B', 'C', 'D'][quad]
		times = times_all[quad]
		lc = lc_all[quad]
		tplot = 0.5 * (times[1:] + times[:-1])
		plt.plot(tplot, lc, alpha=0.8, label=quadname)
		plt.ylabel("Counts/sec")
		plt.xlabel("Time (sec)")
	plt.legend(fontsize='small',handlelength=0.5,labelspacing=0.2,ncols=2)
	plt.title(f"{tbin}s CZT LC")
	for tm in tmark:
		plt.axvline(tm, color='black', alpha=0.8, lw=0.8, linestyle='dashed', label='Marked times')
	plotfile.savefig()
	fig.set_size_inches(5,4)
	fig.set_dpi(120)
	plt.savefig(pdffolder+'/snap_czti_'+str(tmark[0])+'_'+str(tbin)+'_.png')
	plt.close('all')
	plotfile.close()
	old_pdf=PdfReader(plotfilename)
	new_pdf=PdfWriter()
	new_pdf.add_page(old_pdf.pages[-1])
	for i in range(len(old_pdf.pages)-1):
		new_pdf.add_page(old_pdf.pages[i])
	with open(plotfilename,'wb') as f:
		new_pdf.write(f)

	return statfile

def plot_veto_lc(plotflag, statfile, pdffolder, dirname, orbitname, tmark, tbin, tmin, tmax):
	tmark = tmark.astype(float)
	tbin = float(tbin)
	tmin = float(tmin)
	tmax = float(tmax)
	plotfilename = "{stem}_lightcurves_veto.pdf".format(stem=pdffolder+'/'+str(tmark[0])+'_'+str(tbin))
	if(len(glob.glob(plotfilename)) != 0 and plotflag==False):
		return Table.read(pdffolder+'/statfile.csv')
	plotfile = PdfPages(plotfilename)
	eventfile = glob.glob('data/local_level2/'+dirname+'/czti/orbit/'+orbitname+'/*quad_clean.evt')[0]
	#print eventfile
	gains = np.array([5.591,5.594,5.943,5.222])
	offsets = np.array([-56.741, -41.239, -41.682, -26.528])
	channel_emin = 0
	channel_emax = 128
	channel_ebin = 1.0
	hdu = fits.open(eventfile)
	tbl = Table(hdu[5].data)
	tmin = tmin - 50.0*tbin
	tmax = tmax + 50.0*tbin
	channel_ebins = np.arange(channel_emin, channel_emax+channel_ebin/100., channel_ebin)
	num_ebins = len(channel_ebins)-1
	lc_all = []
	times_all = []
	for quad in range(4):
		quadname = ['A', 'B', 'C', 'D'][quad]
		data = tbl[tbl['QuadID']==quad]
		#print data, t_min, t_max, min(data['time']), max(data['time'])
		t_min = max((min(data['Time']), tmin) )  # Start at datastart or args.tmin, whichever comes LATER
		t_max = min((max(data['Time']), tmax) )  # End at dataend or args.tmax, whichever comes FIRST
		hist2d = data['VetoSpec'][:,0:128]
		times = np.arange(t_min, t_max + tbin, tbin)
		times_all.append(times)
		n = int(tbin)
		# print(tmin, min(data['Time']), max(data['Time']), tmark)
		if((len(times) == 0) or tmark[0] >= max(data['Time'])):
			statf_row = [orbitname,'veto',tbin,0,0,0,-100]
			statfile.add_row(statf_row)
			plotfile.close()
			return statfile
		bin_start = np.where(data['Time'] >= t_min)[0][0]
		bin_end = np.where(data['Time'] <= t_max)[0][-1]
		#print bin_start, bin_end
		hist2d = np.array(hist2d[bin_start:bin_end])
		#print len(hist2d)
		ebins = channel_ebins*gains[quad] + offsets[quad]
		emin = channel_emin*gains[quad] + offsets[quad]
		emax = channel_emax*gains[quad] + offsets[quad]
		ebin = int(np.round(channel_ebin*gains[quad],0))
		h = []
		for i in range(len(times)):
			h.append(np.sum(hist2d[n*i:n*(i+1)],axis=0))
			#times.append(np.mean(times_init[i:i+1]))
		h = np.array(h)

		# print(t_min, t_max, tbin)
		# print("veto", np.shape(h), len(times))
		if(len(h)!=len(times)):
			print("length not same error")

		sigmas = np.zeros(num_ebins)
		medians = np.zeros(num_ebins)
		hnorm = np.zeros(np.shape(h))
		hsub = np.zeros(np.shape(h))
		for count in range(num_ebins):
			calc = sigma_clipped_stats(h[:,count], sigma=3.0, maxiters=3)
			medians[count] = calc[1]
			if calc[2] != 0:
				sigmas[count] = calc[2]
			else:
				sigmas[count] = 1.0
			hsub[:,count] = h[:,count] - medians[count]
			hnorm[:,count] = (h[:,count] - medians[count]) / sigmas[count]

		lc = make_plot(plotfile,h,ebin,tbin,quad,emin,emax,ebins,t_min,t_max,times,tmark, filename="{q}_lc_veto".format(q=quadname), title="LC:Quad {q}".format(q=quadname), specscale='log')
		lcsub = make_plot(plotfile,hsub,ebin,tbin,quad,emin,emax,ebins,t_min,t_max,times,tmark, unit="Excess counts", filename="{q}_lcsub_veto".format(q=quadname), title="LC-median subtracted:Quad {q}".format(q=quadname))
		lcnorm = make_plot(plotfile,hnorm,ebin,tbin,quad,emin,emax,ebins,t_min,t_max,times,tmark, unit="Strength", filename="{q}_lcnorm_veto".format(q=quadname), title="LC-normalised:Quad {q}".format(q=quadname))
		#lc_justlast_band = make_plot(plotfile,h[:,110:],ebin,tbin,quad,ebins[110],ebins[-1],ebins[110:],t_min,t_max,times,tmark, filename="{q}_lc_veto".format(q=quadname), title="LC last bands:Quad {q}".format(q=quadname), specscale='log')
		#lcnorm = make_plot(plotfile,hnorm[:,110:],ebin,tbin,quad,ebins[110],ebins[-1],ebins[110:],t_min,t_max,times,tmark, unit="Strength", filename="{q}_lcnorm_veto".format(q=quadname), title="LC-normalised last band:Quad {q}".format(q=quadname))
		lc_all.append(lc)
		
		# Combined lightcurves
		fig = plt.figure()
		gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
		ax_top = fig.add_subplot(gs[0])
		ax_top.set_title("Lightcurve: Quadrant " + quadname)
		tplot = times
		if(len(tplot)!=len(lc)):
			tplot = np.array(0.5 * (tplot[1:] + tplot[:-1]))
		ax_top.plot(tplot, lc, alpha=0.9,label="{emin:0.1f}-{emax:0.1f} keV LC".format(emin=emin, emax=emax))
		ax_top.plot(tplot, lcsub, alpha=0.9,label="Median spectrum subtracted LC")
		ax_top.set_ylabel("Counts/sec", color='black')
		#ax_top.set_xlabel("Time (sec)", color='black')
		ylims = ax_top.get_ylim()
		yupper = ylims[1] + 0.20 * (ylims[1] - ylims[0])
		for tm in tmark:
			ax_top.axvline(tm, color='black', alpha=0.8,lw=0.8,linestyle='dashed', label='Marked times')
		ax_top.set_ylim(ylims[0], yupper)
		ax_top.set_xlim(t_min, t_max)
		plt.xticks(fontsize='small')
		plt.yticks(fontsize='small')
		ax_top.xaxis.set_minor_locator(MultipleLocator(tbin*2.5))	
		ax_top.yaxis.set_minor_locator(AutoMinorLocator(5))
		offset=ax_top.xaxis.get_offset_text()
		offset.set_size('x-small')
		ax_top.legend(fontsize='x-small', handlelength=1, framealpha=0.5)
		ax_top.grid(alpha=0.3)
		#ax_bottom = plt.subplot(gs[1])+
		ax_bottom=fig.add_subplot(gs[1],sharex=ax_top)
		ax_bottom.plot(tplot, lcnorm, label="Mean spectrum subtracted + normalized LC")
		ax_bottom.set_xlabel("Trigger Time (sec)", color='black')
		ax_bottom.set_ylabel("Significance")
		ylims = ax_bottom.get_ylim()
		yupper = ylims[1] + 0.40 * (ylims[1] - ylims[0])
		for tm in tmark:
			ax_bottom.axvline(tm, color='black', alpha=0.8,lw=0.8,linestyle='dashed', label='Marked times')
		ax_bottom.set_ylim(ylims[0], yupper)
		ax_bottom.set_xlim(t_min, t_max)
		plt.xticks(fontsize='small')
		plt.yticks(fontsize='small')
		ax_bottom.xaxis.set_minor_locator(MultipleLocator(tbin*2.5))
		ax_bottom.yaxis.set_minor_locator(AutoMinorLocator(5))
		offset=ax_bottom.xaxis.get_offset_text()
		offset.set_size('x-small')
		ax_bottom.legend(fontsize='x-small', handlelength=1, framealpha=0.5)
		ax_bottom.grid(alpha=0.3)
		plotfile.savefig()
		plt.close()

	t90, bg, rate, counts = calc_stats(lc_all, times_all, tmark, tbin, t_min, t_max, plotfilename, eventfile, plotfile)
	statf_row = [orbitname,'veto',tbin,bg,rate,counts,t90]
	statfile.add_row(statf_row)
	fig = plt.figure()
	for quad in range(4):
		# plot all 4 lightcurves on a single plot
		quadname = ['A', 'B', 'C', 'D'][quad]
		lc = lc_all[quad]
		tplot = np.array(times)
		if(len(tplot)!=len(lc)):
			minlen = min(len(tplot),len(lc))
			if(len(tplot) > len(lc)):
				tplot = tplot[0:minlen]
			else:
				lc = lc[0:minlen]
		plt.plot(tplot, lc, alpha=0.8, label=quadname)
		plt.ylim()
		plt.xlim(t_min,t_max)
		plt.ylabel("Counts/sec")
		plt.xlabel("Time (sec)")
	plt.legend(fontsize='small',handlelength=0.5,labelspacing=0.2,ncols=2)
	plt.title(f"{tbin}s Veto LC")
	for tm in tmark:
		plt.axvline(tm, color='black',alpha=0.8,lw=0.8, linestyle='dashed', label='Marked times')
	plotfile.savefig()
	fig.set_size_inches(5,4)
	fig.set_dpi(120)
	plt.savefig(pdffolder+'/snap_veto_'+str(tmark[0])+'_'+str(tbin)+'_.png')
	plt.close('all')
	plotfile.close()
	old_pdf=PdfReader(plotfilename)
	new_pdf=PdfWriter()
	new_pdf.add_page(old_pdf.pages[-1])
	for i in range(len(old_pdf.pages)-1):
		new_pdf.add_page(old_pdf.pages[i])
	with open(plotfilename,'wb') as f:
		new_pdf.write(f)

	return statfile

def plot_lcs(plotflag, statfile, pdffolder, dirname, orbitname, tmark, all_tbins, tmin, tmax):
	print(f"Plotting for these tbins: {all_tbins}")
	for tbin in all_tbins:
		if(tbin == 0.1):
			print(f'Starting for {tbin} CZT')
			statfile = plot_czti_lc(plotflag, statfile, pdffolder, dirname, orbitname, tmark, tbin, tmin, tmax)
		else:
			print(f'Starting for {tbin} CZT')
			statfile = plot_czti_lc(plotflag, statfile, pdffolder, dirname, orbitname, tmark, tbin, tmin, tmax)
			print(f'Starting for {tbin} Veto')
			statfile = plot_veto_lc(plotflag, statfile, pdffolder, dirname, orbitname, tmark, tbin, tmin, tmax)
	statfile.write(pdffolder+'/statfile.csv', overwrite=True)
	return statfile


# def plot_czti_lc_bandwise(dirname,orbitname,tmark,tbin):
# 	plotfile = PdfPages("{stem}_lightcurves_czti.pdf".format(stem='pdffiles/'+str(tmark)+'/'+str(tmark)+'_'+str(tbin)))
# 	eventfile = glob.glob('data/local_level2/'+dirname+'/czti/orbit/'+orbitname+'/*quad_clean.evt')[0]
# 	emin = [0,50.0,100.0]
# 	emax = [50.0,100.0,200.0]
# 	ebin = [5,5,10]
# 	hdu = fits.open(eventfile)
# 	t_min = tmark - 50.0*tbin
# 	t_max = tmark + 50.0*tbin
# 	for band in range(3):
# 		bandname = ['0-50 keV','50-100 keV','>100 keV'][band]
# 		elower = emin[band]
# 		eupper = emax[band]
# 		e_bin = ebin[band]
# 		ebins = np.arange(elower, eupper+e_bin/100., e_bin)
# 		num_ebins = len(ebins)-1
# 		lc_all = []
# 		times_all = []
# 		for quad in range(4):
# 			quadname = ['A', 'B', 'C', 'D'][quad]
# 			data = hdu[quad+1].data
# 			tmin = max((min(data['time']), t_min) )  # Start at datastart or args.tmin, whichever comes LATER
# 			tmax = min((max(data['time']), t_max) )  # End at dataend or args.tmax, whichever comes FIRST
# 			times = np.arange(tmin, tmax+tbin/100., tbin)
# 			times_all.append(times)
# 			h = np.histogram2d(data['time'], data['energy'], bins=(times, ebins))
# 			h = h[0]
# 			sigmas = np.zeros(num_ebins)
# 			medians = np.zeros(num_ebins)
# 			hnorm = np.zeros(np.shape(h))
# 			hsub = np.zeros(np.shape(h))
# 			for count in range(num_ebins):
# 				calc = sigma_clipped_stats(h[:,count], sigma=3.0, iters=3)
# 				medians[count] = calc[1]
# 				if calc[2] != 0:
# 					sigmas[count] = calc[2]
# 				else:
# 					sigmas[count] = 1.0
# 				hsub[:,count] = h[:,count] - medians[count]
# 				hnorm[:,count] = (h[:,count] - medians[count]) / sigmas[count]

# 			lc = make_plot(plotfile,h,e_bin,tbin,quad,elower,eupper,ebins,tmin,tmax,times,tmark, filename="{q}_lc".format(q=quadname), title="LC:Band {b} Quad {q}".format(b=bandname,q=quadname), specscale='log')
# 			lcsub = make_plot(plotfile,hsub,e_bin,tbin,quad,elower,eupper,ebins,tmin,tmax,times,tmark, unit="Excess counts", filename="{q}_lcsub".format(q=quadname), title="LC-median subtracted:Band {b} Quad {q}".format(b=bandname,q=quadname))
# 			lcnorm = make_plot(plotfile,hnorm,e_bin,tbin,quad,elower,eupper,ebins,tmin,tmax,times,tmark, unit="Strength", filename="{q}_lcnorm".format(q=quadname), title="LC-normalised:Band {b} Quad {q}".format(b=bandname,q=quadname))

# 			lc_all.append(lc)

# 			# Combined lightcurves
# 			fig = plt.figure()
# 			gs = gridspec.GridSpec(2, 1, height_ratios=[2,1])
# 			ax_top = plt.subplot(gs[0])
# 			plt.title("LC:Band "+ bandname+" Quad " + quadname)
# 			tplot = 0.5 * (times[1:] + times[:-1])
# 			plt.plot(tplot, lc, label="{emin:0.1f}-{emax:0.1f} keV lightcurve".format(emin=elower, emax=eupper))
# 			plt.plot(tplot, lcsub, label="Median spectrum subtracted lightcurve")
# 			plt.ylabel("Counts/sec", color='red')
# 			plt.xlabel("Time (sec)", color='red')
# 			ylims = plt.ylim()
# 			yupper = ylims[1] + 0.20 * (ylims[1] - ylims[0])
# 			if tmark > 0: plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
# 			plt.ylim(ylims[0], yupper)
# 			plt.xlim(tmin, tmax)
# 			plt.legend(fontsize='small', framealpha=0.5)
# 			ax_bottom = plt.subplot(gs[1])
# 			plt.plot(tplot, lcnorm, label="Mean spectrum subtracted + normalized lightcurve")
# 			plt.xlabel("Time (sec)", color='red')
# 			plt.ylabel("Significance")
# 			ylims = plt.ylim()
# 			yupper = ylims[1] + 0.40 * (ylims[1] - ylims[0])
# 			if tmark > 0: plt.axvline([tmark], color='red', linestyle='dashed', label='Marked times')
# 			plt.ylim(ylims[0], yupper)
# 			plt.xlim(tmin, tmax)
# 			plt.legend(fontsize='small', framealpha=0.5)
# 			plotfile.savefig()
# 			plt.close()

# 		fig = plt.figure()
# 		for quad in range(4):
# 			# plot all 4 lightcurves on a single plot
# 			quadname = ['A', 'B', 'C', 'D'][quad]
# 			times = times_all[quad]
# 			lc = lc_all[quad]
# 			tplot = 0.5 * (times[1:] + times[:-1])
# 			plt.plot(tplot, lc, label=quadname)
# 			plt.ylabel("Counts/sec")
# 			plt.xlabel("Time (sec)")
# 		plt.legend()
# 		plt.title("All LC band "+bandname)
# 		if tmark > 0: plt.axvline([tmark], color='gray', linestyle='dashed', label='Marked times')
# 		plotfile.savefig()
# 		plt.savefig('pdffiles/'+str(tmark)+'/'+str(tmark)+'_'+str(tbin)+'_czti_band'+str(band)+'.png')
# 		plt.close()
# 	try:
# 		fig = plt.figure()
# 		for quad in range(4):
# 			quadname = ['A', 'B', 'C', 'D'][quad]
# 			veto_file = Table.read('data/local_level2/'+dirname+'/czti/orbit/'+orbitname+'/'+orbitname+'_veto_'+str(tbin)+'_Q'+str(quad)+'_detrended.fits')
# 			times = np.array(veto_file['time'])
# 			tplot = 0.5 * (times[1:] + times[:-1])
# 			tplot = np.concatenate([tplot,[tplot[-1]+tbin]])
# 			#ind_tmin = np.where(abs(tplot - t_min) <= tbin)[0][0]
# 			#ind_tmax = np.where(abs(tplot - t_max) <= tbin)[0][0]
# 			vetolc = np.array(veto_file['vetocounts'])
# 			plt.plot(tplot, vetolc, label = quadname)
# 			plt.xlim(t_min,t_max)
# 		plt.title("Veto lightcurve")
# 		plt.ylabel("Veto Counts/sec")
# 		plt.xlabel("Time (sec)")
# 		plt.legend()	
# 		if tmark > 0: plt.axvline([tmark], color='gray', linestyle='dashed', label='Marked times')
# 		plotfile.savefig()
# 		plt.savefig('pdffiles/'+str(tmark)+'/'+str(tmark)+'_'+str(tbin)+'_veto.png')
# 		plt.close()
# 	except:
# 		print "Error in plot czti"
# 	plotfile.close()
