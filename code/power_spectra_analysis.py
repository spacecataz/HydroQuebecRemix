

import datetime as dt
import numpy as np
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

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__


def open_output_data(outputdir, keywrd, starttime):
    files = sorted(glob.glob(outputdir + '/{}/mag*mag'.format(keywrd)))
    for fpath in files:
        origdata = spacepy.pybats.bats.MagFile(fpath)
        origdata.calc_h()
        origdata.calc_dbdt() 
    return origdata

def build_spectra_plot(ax,fig,date,index):
    ax.legend()
    plt.ylabel(r'$[nT s]^2$')
    ax.set_xlabel('[Hz]')

    ax.axvline(x=0.001111,linestyle='--', color='black')
    ax.axvline(x=0.0005556,linestyle='--', color='black')
    ax.axvline(x=0.0002778,linestyle='--', color='black')

    ax.annotate("15min", xy=[0.001111, 1e4], fontsize=10,rotation=90)
    ax.annotate("30min", xy=[0.0005556, 1e4], fontsize=10,rotation=90)
    ax.annotate("60min", xy=[0.0002778, 1e4], fontsize=10,rotation=90)

    title = 'Power Spectra for {0} for {1}'.format(date, index)
    ax.set_title(title)
    fig.savefig('plots/powerspectra/{0}_combined_{1}_average_powerspectra.png'.format(date,index))


def power_spectra(data,index,date):
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
                print(curdata[:10])
                print(power_spectrum[idx][:10])
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
    file = util.load_logs('/Users/smg3652/Desktop/SWMF_analysis/outputs/{}/{}'.format(date,keywrd),logtype)
    return file

def spectra_per_index(date,index,datalist):
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

    alldata = [['All', 'Obs', obs_data]]
    for keywrd,data in datalist:
        datai = data[index]
        print(index,len(obs_data), len(datai))
        print(data.keys())
        alldata = alldata + [[index,keywrd,datai]]
    power_spectra(alldata,index,date)

def open_all_logs(date,logtype):
    datalist = []
    for keywrd in ['1min','15min','30min','60min']:
        curdata = openlogfile(date,keywrd,logtype)
        datalist = datalist + [[keywrd,curdata]]
    return datalist

starttimes = ['20061214120000','20010831000000','20050831100000','20100405000000','20110805090000']
fig1, ax1 = plt.subplots(figsize=(10,8))

for starttime in starttimes:
    date = starttime[:8]
    if date == '20010831':
        date = '20010830'
    geodatalist=open_all_logs(date,'geo')
    logdatalist=open_all_logs(date,'log')
    if date == '20010830':
        date = '20010831'

    spectra_per_index(date,'AE',geodatalist)
    spectra_per_index(date,'AU',geodatalist)
    spectra_per_index(date,'AL',geodatalist)

    spectra_per_index(date,'dst_sm',logdatalist)

    



