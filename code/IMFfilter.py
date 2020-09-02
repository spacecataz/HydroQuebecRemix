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

def interp_imf(IMF_data):
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

	IMF_data['time'] = newtime
	return IMF_data

def imf_ds(IMF_data, rate):
	# Downsamples 5 second cadence into 5min, 30min, and hourly values
	# scipy decimate does not work well for rates > 13 so we do this progressivly
	print('Beginning downsample')
	downsample = partial(signal.decimate, ftype='iir', axis=-1, zero_phase=True)

	# Downsample to 30 second cadence
	for key in list(IMF_data)[1:]:
		IMF_data[key] = downsample(IMF_data[key],6)
	IMF_data['time'] = IMF_data['time'][::6]

	# Downsample to every 5 min
	for key in list(IMF_data)[1:]:
		IMF_data[key] = downsample(IMF_data[key],10)
	IMF_data['time'] = IMF_data['time'][::10]
	#print(IMF_data['time'][:5])
	IMF_data.quicklook().savefig('new_imf_5min.png')

	if rate == 5:
		return IMF_data

	# Downsample to every 30 min
	for key in list(IMF_data)[1:]:
		IMF_data[key] = downsample(IMF_data[key],6)
	IMF_data['time'] = IMF_data['time'][::6]
	#print(IMF_data['time'][:5])
	IMF_data.quicklook().savefig('new_imf_30min.png')

	if rate == 30:
		return IMF_data
	# Downsample to every hour
	for key in list(IMF_data)[1:]:
		IMF_data[key] = downsample(IMF_data[key],2)
	IMF_data['time'] = IMF_data['time'][::2]
	#print(IMF_data['time'][:5])
	IMF_data.quicklook().savefig('new_imf_hour.png')


	return IMF_data

if __name__ == "__main__":
	rate = int(sys.argv[1])
	fpath = sys.argv[2]
	imf = pybats.ImfInput(fpath)
	new_imf = pybats.ImfInput(fpath)
	interp_imf = interp_imf(new_imf)

	interp_imf.quicklook().savefig('interpolated_imf.png')
	new_imf = imf_ds(interp_imf, rate)
	imf.quicklook().savefig('original_imf.png')