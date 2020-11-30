import datetime as dt
import numpy as np
from scipy import linalg
from scipy.signal import decimate
import matplotlib.pyplot as plt
import spacepy.toolbox as tb
import spacepy.plot as splot
import spacepy.pybats.bats
import sys
sys.path.append('/Users/sgraf/Desktop/SWMFtools')
import util
sys.path.append('/Users/sgraf/Desktop/SWMFtools/dBdt')
import supermag_parser
from spacepy.plot import applySmartTimeTicks
import spacepy.plot as splot
import verify
import os
import glob
splot.style('spacepy')
plt.rcParams["legend.frameon"] = True
plt.rcParams["legend.facecolor"] = 'white'
def read_uncert_line(line, origval, origvalmin, origvalmax):
    tmpline = line.split(' ')
    val = float(tmpline[1])
    valmin = float(tmpline[2][1:-1])
    valmax = float(tmpline[3][:-2])

    origval += [val]
    origvalmin += [val-valmin]
    origvalmax += [valmax-val]

    return origval, origvalmin, origvalmax

def read_file(fname, thresholds):
    with open(fname,'r') as f:
        lines = f.readlines()
        POFD, POFDmin, POFDmax = [], [], []
        POD, PODmin, PODmax = [], [], []
        h, hmin, hmax = [], [], []
        bias, biasmin, biasmax = [], [], []

        for threshold in thresholds:
            while lines[0][0] != 'T': lines = lines[1:]
            threshold = float(lines[0].split(':')[1][:-1])
            while lines[0][:4] != 'Bias': lines = lines[1:]
            bias, biasmin, biasmax = read_uncert_line(lines[0], bias, biasmin, biasmax)
            while lines[0][:4] != 'POFD': lines = lines[1:]
            POFD, POFDmin, POFDmax = read_uncert_line(lines[0], POFD, POFDmin, POFDmax)
            POD, PODmin, PODmax = read_uncert_line(lines[1], POD, PODmin, PODmax)
            while lines[0][:4] != 'Heid': lines = lines[1:]
            h, hmin, hmax = read_uncert_line(lines[0], h, hmin, hmax)

        POFDerr = np.zeros((2,4))
        POFDerr[0,] = np.array(POFDmin)
        POFDerr[1,] = np.array(POFDmax)

        PODerr = np.zeros((2,4))
        PODerr[0,] = np.array(PODmin)
        PODerr[1,] = np.array(PODmax)

        herr = np.zeros((2,4))
        herr[0,] = np.array(hmin)
        herr[1,] = np.array(hmax)

        biaserr = np.zeros((2,4))
        biaserr[0,] = np.array(biasmin)
        biaserr[1,] = np.array(biasmax)
    f.close()
    return POFD, POFDerr, POD, PODerr, h, herr, bias, biaserr
def main(args):
    date = args[1]
    unsmoothed_fnames = glob.glob('ctables/{}*unsmoothed.txt'.format(date))

    thresholds = [0.3, 0.7, 1.1, 1.5] #nT/s #nT/s
    thresholds_unsmoothed = [0.28, 0.68, 1.08, 1.48]
    thresholds_30min = [0.3, 0.7, 1.1, 1.5] #nT/s
    thresholds_hour = [0.32, 0.72, 1.12, 1.52] #nT/s
    for fname in unsmoothed_fnames:
        hourly_fname = fname[:-14] + 'hour.txt'
        fname_30 = fname[:-14] + '30min.txt'
        sat = fname.split('/')[-1].split('_')[1]
        POFD, POFDerr, POD, PODerr, h, herr, bias, biaserr = read_file(fname, thresholds)
        hour_POFD, hour_POFDerr, hour_POD, hour_PODerr, hour_h, hour_herr, hour_bias, hour_biaserr = read_file(hourly_fname, thresholds)
        POFD_30, POFDerr_30, POD_30, PODerr_30, h_30, herr_30, bias_30, biaserr_30 = read_file(fname_30, thresholds)
        
        '''
        # pulk plot type
        fig, axs = plt.subplots(figsize=(10,6),nrows=4, ncols=2, sharex=True)
        for i in range(4):
            axs[i,0].errorbar(1, POFD[i], yerr=POFDerr[:,i], fmt='bo',  label = 'Unsmoothed')
            axs[i,0].errorbar(1, POD[i], fmt='ko', yerr = PODerr[i], label = 'Unsmoothed')
            axs[i,1].errorbar(1, h[i], fmt='o', yerr = herr[i], label = 'Unsmoothed')   

            axs[i,0].errorbar(2, hour_POFD[i], fmt='bo', yerr = hour_POFDerr[i], label = 'Unsmoothed')
            axs[i,0].errorbar(2, hour_POD[i], fmt='ko', yerr = hour_PODerr[i], label = 'Unsmoothed')
            axs[i,1].errorbar(2, hour_h[i], fmt='o', yerr = hour_herr[i], label = 'Unsmoothed')      

        plt.show()
        '''
        # plot type 1
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(figsize=(10,12),nrows=4, sharex=True)

        print(sat, h)
        ax1.errorbar(thresholds_unsmoothed, POFD, fmt='o', yerr = POFDerr, capsize=3, label = 'Unsmoothed')
        ax2.errorbar(thresholds_unsmoothed, POD, fmt='o', yerr = PODerr, capsize=3, label = 'Unsmoothed')
        ax3.errorbar(thresholds_unsmoothed, h, fmt='o', yerr = herr, capsize=3, label = 'Unsmoothed')
        ax4.errorbar(thresholds_unsmoothed, bias, fmt='o', yerr = biaserr, capsize=3, label = 'Unsmoothed')

        ax1.errorbar(thresholds_30min, POFD_30, fmt='s', yerr = POFDerr_30, capsize=3, label = '30min')
        ax2.errorbar(thresholds_30min, POD_30, fmt='s', yerr = PODerr_30, capsize=3, label = '30min')
        ax3.errorbar(thresholds_30min, h_30, fmt='s', yerr = herr_30, capsize=3, label = '30min')
        ax4.errorbar(thresholds_30min, bias_30, fmt='s', yerr = biaserr_30, capsize=3, label = '30min')

        ax1.errorbar(thresholds_hour, hour_POFD, fmt='s', yerr = hour_POFDerr, capsize=3, label = 'Hourly')
        ax2.errorbar(thresholds_hour, hour_POD, fmt='s', yerr = hour_PODerr, capsize=3, label = 'Hourly')
        ax3.errorbar(thresholds_hour, hour_h, fmt='s', yerr = hour_herr, capsize=3, label = 'Hourly')
        ax4.errorbar(thresholds_hour, hour_bias, fmt='s', yerr = hour_biaserr, capsize=3, label = 'Hourly')




       
        ax1.set_title('POFD')
        ax2.set_title('POD')
        ax3.set_title('Heidke Score')
        ax4.set_title('Bias')

        ax1.legend()
        #ax2.legend()
        #ax3.legend()
        #ax4.legend()

        ax1.set_ylim([0,1])
        ax2.set_ylim([0,1])
        ax3.set_ylim([-0.2,1])

        ax4.set_xlim([0,1.8])

        ax1.set_xticks(thresholds)
        ax4.set_xlabel('Thresholds [nT/s]')
        ax3.plot([0,1.8], [0,0], 'k--')
        ax4.plot([0,1.8], [1,1], 'k--')

        fig.suptitle('Model Statistics Summary for {}'.format(sat))
        fig.tight_layout()
        fig.subplots_adjust(top=0.92)
        plt.savefig('plots/{0}_{1}_skill_comp.png'.format(date,sat))
        #plt.show()




if __name__ == "__main__":
    main(sys.argv)
