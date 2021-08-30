#run by python IMFfilter.py [downsample rate in minutes] [path to file]
#input oversampled to 5 second cadence

import spacepy
import spacepy.pybats as pybats
from functools import partial
from scipy import signal
import sys
from scipy.interpolate import interp1d
import datetime
import matplotlib.pyplot as plt
from datetime import datetime as dt
import glob
import spacepy.plot as splot
splot.style('spacepy')

def interp_imf(IMF_data, path):
	# use this to create evenly spaced time steps, we will oversample here
	newtime = []
	for key in list(IMF_data)[1:]:
		print('Interpolating for {}'.format(key))
		timestamps = [y.timestamp() for y in IMF_data['time']]
		k1 = interp1d(timestamps,IMF_data[key],'linear')
		starttime = IMF_data['time'][0]
		endtime = IMF_data['time'][-1]
		startlist = [starttime.replace(microsecond=0) + datetime.timedelta(seconds=1)]


		while startlist[-1] < endtime - datetime.timedelta(seconds=5):
			startlist = startlist + [startlist[-1] + datetime.timedelta(seconds=5)]

		eventimestamps = [t.timestamp() for t in startlist]

		new_imf_data = k1(eventimestamps)
		IMF_data[key] = new_imf_data
		newtime = startlist
	print('Writing IMF')
	IMF_data['time'] = newtime
	IMF_data.write(outfile=path + 'IMF_interpolated.dat')
	return IMF_data

def imf_ds(IMF_data, rate, path):
	# Downsamples 5 second cadence into 5min, 30min, and hourly values
	# scipy decimate does not work well for rates > 13 so we do this progressivly

	# try subtracting the mean




	print('Beginning downsample')
	downsample = partial(signal.decimate, ftype='iir', axis=-1, zero_phase=True)

	# Downsample to 30 second cadence
	for key in list(IMF_data)[1:]:
		mean = IMF_data[key].mean()
		data_sub = IMF_data[key] - mean
		down_sub = downsample(data_sub,6)
		print(mean, down_sub[0])
		IMF_data[key] = down_sub + mean
		print(IMF_data[key][0])
		print('HI')

	#times = [dt.fromtimestamp(t) for t in newtimestamps]
	times = IMF_data['time'][::6]

	IMF_data['time'] = times

	# Downsample to every 5 min
	for key in list(IMF_data)[1:]:
		#IMF_data[key] = downsample(IMF_data[key],10)
		mean = IMF_data[key].mean()
		data_sub = IMF_data[key] - mean
		down_sub = downsample(data_sub,10)
		IMF_data[key] = down_sub + mean

	times = IMF_data['time'][::10]

	IMF_data['time'] = times


	IMF_data.write(outfile=path + 'IMF_5min.dat')

	if rate == 5:
		return IMF_data

	for key in list(IMF_data)[1:]:
		#IMF_data[key] = downsample(IMF_data[key],6)
		mean = IMF_data[key].mean()
		data_sub = IMF_data[key] - mean
		down_sub = downsample(data_sub,3)
		IMF_data[key] = down_sub + mean

	times = IMF_data['time'][::3]

	IMF_data['time'] = times


	IMF_data.write(outfile=path + 'IMF_15min.dat')

	# Downsample to every 30 min
	for key in list(IMF_data)[1:]:
		#IMF_data[key] = downsample(IMF_data[key],6)
		mean = IMF_data[key].mean()
		data_sub = IMF_data[key] - mean
		down_sub = downsample(data_sub,3)
		IMF_data[key] = down_sub + mean

	times = IMF_data['time'][::3]

	IMF_data['time'] = times


	IMF_data.write(outfile=path + 'IMF_30min.dat')

	if rate == 30:
		return IMF_data
	# Downsample to every hour
	for key in list(IMF_data)[1:]:
		#IMF_data[key] = downsample(IMF_data[key],2)
		mean = IMF_data[key].mean()
		data_sub = IMF_data[key] - mean
		down_sub = downsample(data_sub,2)
		IMF_data[key] = down_sub + mean


	times = IMF_data['time'][::2]

	IMF_data['time'] = times


	IMF_data.write(outfile=path + 'IMF_hourly.dat')
	#print(IMF_data['time'][:5])


	return IMF_data

def plot_for_dan(IMF_data_list):
	import matplotlib.pyplot as plt
	from spacepy.plot import applySmartTimeTicks
	        
	timerange = [dt(2006,12,14,4,12,00), dt(2006,12,16,00,00,00)]

	def adjust_plots(ax, ylab, xlab=False, Zero=True):
	    ax.grid(True)
	    applySmartTimeTicks(ax,timerange)
	    ax.set_ylabel(ylab)
	    labels =ax.get_yticklabels()
	    labels[-1].set_visible(False)
	    labels[0].set_visible(False)
	    if Zero:
	         ax.plot(timerange, [0,0], 'k--')
	    if xlab:
	         ax.set_xlabel('Universal Time from %s' % 
	                        timerange[0].isoformat())
	    else:
	         ax.xaxis.set_ticklabels([])

	fig = plt.figure(figsize=(8,8))
	fig.subplots_adjust(hspace=0.025, top=0.95, bottom=0.05, right=0.95)

	a3 = fig.add_subplot(311)
	a3.plot(IMF_data_list[0]['time'], IMF_data_list[0]['bz'], lw=1.25, c='steelblue', label='Interpolated')
	#a3.plot(IMF_data_list[1]['time'], IMF_data_list[1]['bz'], lw=1.25, c='steelblue', label='5min')
	a3.plot(IMF_data_list[2]['time'], IMF_data_list[2]['bz'], lw=1.25, c='navy', label='30min')
	a3.plot(IMF_data_list[3]['time'], IMF_data_list[3]['bz'], lw=1.25, c='orange',label='Hour')
	#leg = a3.legend()

	adjust_plots(a3, 'IMF $B_{Z}$ ($nT$)')

	a4 = fig.add_subplot(312)
	a4.plot(IMF_data_list[0]['time'], IMF_data_list[0]['rho'], lw=1.25, c='steelblue',label='Interpolated')
	#a4.plot(IMF_data_list[1]['time'], IMF_data_list[1]['rho'], lw=1.25, c='steelblue',label='5min')
	a4.plot(IMF_data_list[2]['time'], IMF_data_list[2]['rho'], lw=1.25, c='navy', label='30min')
	a4.plot(IMF_data_list[3]['time'], IMF_data_list[3]['rho'], lw=1.25, c='orange', label='hour')
	leg = a4.legend()

	adjust_plots(a4, 'Density ($cm^{-3}$)', Zero=False)

	a5 = fig.add_subplot(313)
	a5.plot(IMF_data_list[0]['time'], -1.0*IMF_data_list[0]['ux'], lw=1.25, c='steelblue', label='Interpolated')
	#a5.plot(IMF_data_list[1]['time'], -1.0*IMF_data_list[1]['ux'], lw=1.25, c='steelblue', label='5min')
	a5.plot(IMF_data_list[2]['time'], -1.0*IMF_data_list[2]['ux'], lw=1.25, c='navy', label='30min')
	a5.plot(IMF_data_list[3]['time'], -1.0*IMF_data_list[3]['ux'], lw=1.25, c='orange', label='hour')
	#leg = a5.legend()


	plt.tight_layout()
	adjust_plots(a5, '$|V_{X}|$ ($km/s$)', Zero=False, xlab=True)

	return fig

def process_imf(fpath, rate):
	orig_imf = pybats.ImfInput(fpath)
	interped_imf = interp_imf(orig_imf)
	proc_imf = 	imf_ds(interped_imf, rate)

	return proc_imf

if __name__ == "__main__":
	rate = int(sys.argv[1])
	fpath = sys.argv[2]
	path = fpath[:-7]
	print(path)
	imf = pybats.ImfInput(fpath,path)
	new_imf = pybats.ImfInput(fpath,path)
	interped_imf = interp_imf(new_imf,path)
	new_imf = imf_ds(interped_imf, rate,path)

	five_imf = 	pybats.ImfInput('/Users/sgraf/Desktop/HydroQuebecRemix/swmf_input/Event20061214/IMF_5min.dat')
	thirty_imf = 	pybats.ImfInput('/Users/sgraf/Desktop/HydroQuebecRemix/swmf_input/Event20061214/IMF_30min.dat')
	hour_imf = pybats.ImfInput('/Users/sgraf/Desktop/HydroQuebecRemix/swmf_input/Event20061214/IMF_hourly.dat')
	interped_imf = pybats.ImfInput('/Users/sgraf/Desktop/HydroQuebecRemix/swmf_input/Event20061214/IMF_interpolated.dat')

	fig = plot_for_dan([interped_imf,five_imf,thirty_imf,hour_imf])
	
	plt.savefig('Event6_summary.png')

	#imf.quicklook().savefig('original_imf.png')



