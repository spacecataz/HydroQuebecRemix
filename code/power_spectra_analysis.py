
import spacepy
import spacepy.pybats as pybats
import datetime as dt
import numpy as np
import scipy
from scipy import linalg
from scipy.signal import decimate
import matplotlib.pyplot as plt
import spacepy.toolbox as tb
import spacepy.plot as splot
import spacepy.pybats.bats
import sys
sys.path.append('/Users/smg3652/python_packages/SWMFtools')
import util
sys.path.append('/Users/smg3652/python_packages/SWMFtools/dBdt')
import supermag_parser
from spacepy.plot import applySmartTimeTicks
import verify
import os
import matplotlib.colors as colors
import glob
import matplotlib
from scipy.interpolate import interp1d
import datetime

# This script calculates power spectra plots for inputs + outputs, other than dBdt which is handled by validate_script.py

def open_output_data(outputdir, keywrd, starttime):
    files = sorted(glob.glob(outputdir + '/{}/mag*mag'.format(keywrd)))
    for fpath in files:
        origdata = spacepy.pybats.bats.MagFile(fpath)
        origdata.calc_h()
        origdata.calc_dbdt() 
    return origdata

def process_dBdth(origdata, stat):
    # calculates horizontal dBdt
    dBdth_total = []
    simtime_total = []
    subset = origdata[stat]
    simtime = subset['time'][::6]
    dBdte = decimate(subset['dBdte'], 6)
    dBdtn = decimate(subset['dBdtn'], 6)
    dBdth = np.array([linalg.norm([dBdtn[i], dBdte[i]]) for i in range(len(dBdtn))])
    #dBdth = np.array([linalg.norm([dBdtn[i], dBdte[i]], axis = 0)])
    if len(simtime_total) == 0:
        dBdth_total = dBdth
        simtime_total = simtime

    else:
        dBdth_total = np.concatenate((dBdth_total, dBdth))
        simtime_total = np.concatenate((simtime_total, simtime))
    return dBdth_total

def process_dBh(origdata, stat):
    # calculates horizontal B
    subset = origdata[stat]
    dBe = decimate(subset['dBe'], 6)
    dBn = decimate(subset['dBn'], 6)
    Bh = np.array([linalg.norm([dBn[i], dBe[i]]) for i in range(len(dBn))])
    datamin, datatmin = tb.windowMean(Bh, time=subset['time'][::6], winsize=dt.timedelta(minutes=1), overlap=dt.timedelta(0), st_time=starttime, op=np.max)

    return datamin=

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

def build_spectra_plot(ax,fig,date,index,window_type='None'):
    '''
    Builds formatting of spectra plot.

    INPUTS
    ax, fig: ax, fig of the spectra plot
    date: string date of event ('20061214')
    index: index that the power spectra is of ('bz', 'AU')

    OUTPUTS
    saves resulting figure
    '''

    ax.legend()

    # units
    plt.ylabel(r'$[nT s]^2$')
    ax.set_xlabel('[Hz]')

    # add in 15 min, 30 min, 60 min markers
    ax.axvline(x=0.001111,linestyle='--', color='black')
    ax.axvline(x=0.0005556,linestyle='--', color='black')
    ax.axvline(x=0.0002778,linestyle='--', color='black')

    ax.annotate("15min", xy=[0.001111, 1e4], fontsize=10,rotation=90)
    ax.annotate("30min", xy=[0.0005556, 1e4], fontsize=10,rotation=90)
    ax.annotate("60min", xy=[0.0002778, 1e4], fontsize=10,rotation=90)

    title = 'Power Spectra for {0} for {1}'.format(date, index)
    ax.set_title(title)

    if window_type != 'None':
        fig.savefig('plots/powerspectra/{0}_combined_{1}_average_{2}_powerspectra.png'.format(date,index,window_type))
    else:
        fig.savefig('plots/powerspectra/{0}_combined_{1}_average_powerspectra.png'.format(date,index))


def welch_spectra(data,index,date,window_type):
    '''
    Calculates welch powerspectra of data. 

    INPUTS:
    data: datalist of form ['all', smoothing level ('1min'), data]
    index: index that the power spectra is being calculated of
    date: string date of event ('20061214')
    window_type: type of window for welch processing ('flattop', 'boxcar')
    '''

    fig, ax = plt.subplots(figsize=(10,8))
    return_power = []
    return_freq = []
    if len(data)>3: # this is the toggle to plot all smoothing levels + obs on one plot
        for i in range(len(data)): # loop through out list of data
            # in this case each data entry is: [all/sector, obs or smoothing level, data]
            
            curdata = np.array(data[i][2]) # data
            keywrd = data[i][1] # smoothing level 

            # calculate welch spectra
            frequency, power_spectrum = (scipy.signal.welch(curdata, fs=1.0/60,
                                                            window=scipy.signal.get_window(window_type, int(len(curdata)/5)),scaling='density'))
            
            # cut off first point, because it was causing scales to be massively off...
            power_spectrum = power_spectrum[1:]
            frequency = frequency[1:]

            if i ==0: # plot observed
                plt.loglog(frequency, power_spectrum,'k',label = keywrd)
            else: # plot outputs
                plt.loglog(frequency, power_spectrum,label = keywrd)
            return_power += [np.array(power_spectrum)]
            return_freq += [np.array(frequency)]

        # format plot before saving
        build_spectra_plot(ax,fig,date,index,'welch_{}'.format(window_type))

    else: # this is the case where we are plotting just one spectra
        curdata = np.array(data[2])

        # calculate welch spectra
        frequency, power_spectrum = (scipy.signal.welch(curdata, fs=1.0/60,
                                                        window=scipy.signal.get_window(window_type, int(len(curdata)/5)),scaling='density'))
        
        # cut off first point, because it was causing scales to be massively off...
        power_spectrum = power_spectrum[1:]
        frequency = frequency[1:]

        plt.loglog(frequency, power_spectrum)

        return_power = [power_spectrum[idx]]
        return_freq = [frequency[idx]]

        build_spectra_plot(ax,fig,date,index,'welch_{}'.format(window_type))

    #plt.show()
    fig.savefig('plots/powerspectra/{}/{}_{}_{}_welch_{}_powerspectra.png'.format(date,date,index,data[0][0],window_type))
    plt.close()

    #return return_freq, return_power


def power_spectra(data,index,date):
    '''
    Calculates rfft powerspectra of data. 

    INPUTS:
    data: datalist of form ['all', smoothing level ('1min'), data]
    index: index that the power spectra is being calculated of
    date: string date of event ('20061214')

    '''
    fig, ax = plt.subplots(figsize=(10,8))
    return_power = []
    return_freq = []
    if len(data)==5:
        for i in range(5):
            curdata = np.array(data[i][2])
            keywrd = data[i][1]
            curdata = np.array(curdata)
            
            fourier_transform = np.fft.rfft(curdata)
            abs_fourier_transform = np.abs(fourier_transform)
            power_spectrum = np.square(abs_fourier_transform)
            time_step = 60
            frequency = np.fft.rfftfreq(curdata.size, time_step)
            idx = np.argsort(frequency)

            if i ==0:
                plt.loglog(frequency[idx], power_spectrum[idx],'k',label = keywrd)
            else:
                plt.loglog(frequency[idx], power_spectrum[idx],label = keywrd)
            return_power += [np.array(power_spectrum[idx])]
            return_freq += [np.array(frequency[idx])]


        build_spectra_plot(ax,fig,date,index)
        plt.legend()
    else:
        data = np.array(data[2])
        fourier_transform = np.fft.rfft(data)
        abs_fourier_transform = np.abs(fourier_transform)
        power_spectrum = np.square(abs_fourier_transform)
        time_step = 60
        frequency = np.fft.rfftfreq(data.size, time_step)
        idx = np.argsort(frequency)

        plt.loglog(frequency[idx], power_spectrum[idx])

        return_power = [power_spectrum[idx]]
        return_freq = [frequency[idx]]


    plt.ylabel(r'$[nT s]^2$')
    plt.xlabel('[Hz]')

    #plt.show()

    fig.savefig('plots/powerspectra/{}/{}_{}_{}_powerspectra.png'.format(date,date,index,data[0][0]))
    plt.close()
    return return_freq, return_power

def openlogfile(date,keywrd,logtype):
    '''
    Opens a log type output file

    INPUTS
    date: string date of event ('20061214')
    keywrd: smoothing level ('1min')
    logtype: type of log file ('geo') 

    RETURNS
    file: loaded log instance
    '''
    file = util.load_logs('/Users/smg3652/Desktop/SWMF_analysis/outputs/{}/{}'.format(date,keywrd),logtype)
    return file

def openinputfile(date,keywrd):
    '''
    Opens a IMF input file

    INPUTS
    date: string date of event ('20061214')
    keywrd: smoothing level ('1min')


    RETURNS
    data: IMF data dictionary
    '''
    shortdate = date[:8]
    event = 'Event{}'.format(shortdate)

    # due to naming convention
    if keywrd == '1min':
        keywrd = 'interpolated'

    data = pybats.ImfInput('/Users/smg3652/Desktop/HydroQuebecRemix/swmf_input/{}/IMF_{}.dat'.format(event,keywrd))
    
    return data



def interp_imf(IMF_data):
    '''
    Interpolates data to even timestamps (involves over sampling the data)

    INPUTS
    IMF_data: IMF data dictionary

    RETURNS
    IMF_data: IMF data interpolated to 15 second cadence
    '''

    newtime = [] # list of new timestamps
    for key in list(IMF_data)[1:]:

        timestamps = [y.timestamp() for y in IMF_data['time']] # calculate timestamp for each datetime value
        k1 = interp1d(timestamps,IMF_data[key],'linear') # linearly interpolate between all data points and timestamps
        starttime = IMF_data['time'][0]
        endtime = IMF_data['time'][-1]

        startlist = [starttime.replace(microsecond=0) + datetime.timedelta(seconds=1)] # start to build our final time list
        while startlist[-1] < endtime - datetime.timedelta(seconds=15):
            startlist = startlist + [startlist[-1] + datetime.timedelta(seconds=15)] # add in time steps every 15 seconds

        eventimestamps = [t.timestamp() for t in startlist] # convert to timestamps to index into k1, our linear interpolated data
        new_imf_data = k1(eventimestamps)

        IMF_data[key] = new_imf_data # replace old data spot with new, interpolated data
        newtime = startlist # keep track of new time list, this is the same for every index 

    IMF_data['time'] = newtime # update time entry in IMF_data

    return IMF_data

def spectra_per_input(date,index,datalist, window_type):
    '''
    Function to calculate power spectra for a given input type

    INPUTS
    date: string date of event ('20061214')
    index: index
    datalist: list of data in form [[smoothing level, data],[smoothing level, data],...]
    window_type: type of window for welch procesing ('boxcar', 'flattop') or list of window types

    RETURNS
    nothing, calls welch_spectra which saves the plots
    '''

    # start our data list, in this case we dont have observed so it starts empty
    alldata = []

    # loop through smoothing level, data pairs to build final list
    for keywrd,data in datalist:
        datai = data[index]
        alldata = alldata + [[index,keywrd,datai]]

    # calc power spectra
    if type(window_type)==list:
        for window in window_type:
            welch_spectra(alldata,index,date,window)
    else:
        welch_spectra(alldata,index,date,window_type)

def spectra_per_index(date,index,datalist,window_type):
    '''
    Function to calculate power spectra for a given output index
    INPUTS
    date: string date of event ('20061214')
    index: index
    datalist: list of data in form [[smoothing level, data],[smoothing level, data],...]
    window_type: type of window for welch procesing ('boxcar', 'flattop') or list of window types

    RETURNS
    nothing, calls welch_spectra which saves the plots
    '''

    # grab observed data of the index
    if index == 'AE':
        datalist[0][1].fetch_obs_ae()
        obs_data = datalist[0][1].obs_ae['ae']
    if index == 'AU':
        datalist[0][1].fetch_obs_ae()
        obs_data = datalist[0][1].obs_ae['au']
    if index == 'AL':
        datalist[0][1].fetch_obs_ae()
        obs_data = datalist[0][1].obs_ae['al']
    if index == 'dst_sm':
        datalist[0][1].fetch_obs_sym()
        obs_data = datalist[0][1].obs_sym['sym-h']

    # start our data list
    alldata = [['All', 'Obs', obs_data]] # this goes by type (all or index), obs or smoothing level, data

    # build out our data list for all levels
    for keywrd,data in datalist:
        datai = data[index]

        # this is because of a mismatch in cadence between our ouput and observed dst
        if index == 'dst_sm':
            datai = datai[::12]

        alldata = alldata + [[index,keywrd,datai]]

    # calc power spectra
    if type(window_type)==list:
        for window in window_type:
            welch_spectra(alldata,index,date,window)
    else:
        welch_spectra(alldata,index,date,window_type)
    

def open_all_logs(date,logtype):
    '''
    Opens logfiles for all smoothing levels for a given date

    INPUTS
    date: date of event ('20061214')
    logtype: geo or log

    RETURNS
    datalist: List containing [smoothing level, data] pairs for all levels
    '''
    datalist = [] # this is where we will store our smoothing level, data pairs

    # loop through smoothing levels
    for keywrd in ['1min','15min','30min','60min']:
        curdata = openlogfile(date,keywrd,logtype)
        datalist = datalist + [[keywrd,curdata]]
    return datalist

def open_all_inputs(date):
    '''
    Opens input files for all smoothing levels for a given date

    INPUTS
    date: date of event ('20061214')

    RETURNS
    datalist: List containing [smoothing level, data] pairs for all levels
    '''

    datalist =[] # this is where we will store our smoothing level, data pairs

    # loop through smoothing levels
    for keywrd in ['1min','15min','30min','60min']:

        # this is also because of a naming convention in my SWMF analysis folder
        if keywrd == '60min':
            keywrd = 'hourly'
        curdata = openinputfile(date,keywrd)
        curdata = interp_imf(curdata)
        if keywrd == 'hourly':
            keywrd = '60min'
        datalist = datalist + [[keywrd,curdata]]
    return datalist

# MAIN SECTION OF THE CODE TO RUN IT 

# list of start times to analyze
starttimes = ['20061214120000','20010831000000','20050831100000','20100405000000','20110805090000']
#fig1, ax1 = plt.subplots(figsize=(10,8))

# loop through start times
for starttime in starttimes:
    date = starttime[:8]

    # this is because of a naming convention I have in my SWMF analysis folder...
    if date == '20010831':
        date = '20010830'
    geodatalist=open_all_logs(date,'geo')
    logdatalist=open_all_logs(date,'log')
    if date == '20010830':
        date = '20010831'
    inputdatalist = open_all_inputs(date)

    # calculate power spectra for input indices
    input_indices = ['bz', 'rho', 'ux']
    for index in input_indices:
        spectra_per_input(date,index,inputdatalist,'flattop')

    # calculate power spectra for output indices
    output_indices = ['AE', 'AU', 'AL', 'dst_sm']
    for index in output_indices:
        spectra_per_index(date,index,geodatalist,'flattop')


    



