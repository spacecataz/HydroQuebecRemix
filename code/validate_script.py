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
import verify
import os
import glob

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

def main(args):

    if len(args) >=3:
        date = args[1]
        starttime = args[2]
    else:
        starttime = args[1]
        date = starttime[:8]
    year = date[:4]
    month = date[4:6]
    day = date[6:8]
    hour = starttime[8:10]
    minute = starttime[10:12]
    second = starttime[12:]
    starttime = dt.datetime.strptime(starttime, '%Y%m%d%H%M%S')

    blockPrint()
    smdata = supermag_parser.supermag_parser('./supermag_data/{0}{1}{2}-supermag.txt'.format(year,month,day))
    enablePrint()
    stations = smdata.station
    #stations = [key for key in origdata.keys() if len(key)==3]

    outputdir = './outputs/{}/'.format(date)

    thresholds = [0.3, 0.7, 1.1, 1.5] #nT/s

    #make_table(smdata,stations,outputdir,date,starttime, thresholds,'unsmoothed')
    #make_table(smdata,stations,outputdir,date,starttime, thresholds,'hour')
    #make_table(smdata,stations,outputdir,date,starttime, thresholds,'30min')

    midlat_sats = ['OTT', 'NEW', 'WNG', 'MEA']
    highlat_sats = ['ABK', 'YKC', 'IQA']
    print('high lat')
    #grouping(outputdir, smdata, thresholds, highlat_sats, date, starttime)
    print('midlat')
    #grouping(outputdir, smdata, thresholds, midlat_sats, date, starttime)

    starttimes = ['20061214120000','20010831000000','20050831100000','20100405000000','20110805090000']

    #cross_event_grouping(outputdir, thresholds, highlat_sats, starttimes)
    cross_event_grouping(outputdir, thresholds, midlat_sats, starttimes)
    #grouping(outputdir, smdata, thresholds, midlat_sats, 'Mid_Latitude', date, starttime)
    #grouping(outputdir, smdata, thresholds, highlat_sats, 'High_Latitude', date, starttime)


def grouping(outputdir, smdata, thresholds, stations, date, starttime):
    for threshold in thresholds:
        for keywrd in ['hour', 'unsmoothed', '30min']:
            predicted_event_tot = []
            obs_event_tot = []
            for stat in stations:
                if stat not in smdata.station: continue
                else:
                    dBdth_total, simtime_total = read_output_data(outputdir, keywrd, stat, starttime)
                    smstat = smdata.station[stat]
                    predicted_event_sat, obs_event_sat = make_boolean_arrays(dBdth_total, smstat, simtime_total, starttime, threshold)
                    predicted_event_tot += predicted_event_sat.tolist()
                    obs_event_tot += obs_event_sat.tolist()

            ctable = verify.Contingency2x2.fromBoolean(predicted_event_tot, obs_event_tot)
            if stations[0] == 'ABK':
                group = 'highlat'
            else:
                group = 'midlat'
            write_table(ctable, date, group, keywrd, threshold)

def cross_event_grouping(outputdir, thresholds, stations, starttimes):
    for threshold in thresholds:
        for keywrd in ['hour', 'unsmoothed', '30min']:   
            predicted_event_tot = []
            obs_event_tot = []
            for starttime in starttimes:
                date = starttime[:8]
                starttime = dt.datetime.strptime(starttime, '%Y%m%d%H%M%S')
                if date == '20010831': date = '20010830'
                outputdir = './outputs/{}/'.format(date)
                blockPrint()
                smdata = supermag_parser.supermag_parser('./supermag_data/{0}-supermag.txt'.format(date))
                enablePrint()
                for stat in stations:
                    if stat not in smdata.station: continue
                    else:
                        dBdth_total, simtime_total = read_output_data(outputdir, keywrd, stat, starttime)
                        smstat = smdata.station[stat]
                        predicted_event_sat, obs_event_sat = make_boolean_arrays(dBdth_total, smstat, simtime_total, starttime, threshold)
                        predicted_event_tot += predicted_event_sat.tolist()
                        obs_event_tot += obs_event_sat.tolist()
            print('Lengths before CTable')
            print(len(predicted_event_tot), len(obs_event_tot))
            ctable = verify.Contingency2x2.fromBoolean(predicted_event_tot, obs_event_tot)
            if stations[0] == 'ABK':
                group = 'combined_highlat'
            else:
                group = 'combined_midlat'
            write_table(ctable, date, group, keywrd, threshold)

def read_output_data(outputdir, keywrd, stat, starttime):
    files = sorted(glob.glob(outputdir + '/{}/mag*mag'.format(keywrd)))
    dBdth_total = []
    simtime_total = []
    for fpath in files:
        origdata = spacepy.pybats.bats.MagFile(fpath)
        origdata.calc_h()
        origdata.calc_dbdt()

        subset = origdata[stat]
        simtime = subset['time'][::6]
        dBdte = decimate(subset['dBdte'], 6)
        dBdtn = decimate(subset['dBdtn'], 6)
        dBdth = np.array([linalg.norm([dBdtn[i], dBdte[i]]) for i in range(len(dBdtn))])

        if len(simtime_total) == 0:
            dBdth_total = dBdth
            simtime_total = simtime

        else:
            dBdth_total = np.concatenate((dBdth_total, dBdth))
            simtime_total = np.concatenate((simtime_total, simtime))
    return dBdth_total, simtime_total

def make_boolean_arrays(dBdth, smstat, simtime, starttime, threshold):
    run20, t20 = tb.windowMean(dBdth, time=simtime, winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
    predicted_event = np.asarray(run20) >= threshold

    Bdoth = np.array([linalg.norm(smstat['Bdot'][i,:2]) for i in range(len(smstat['Bdot']))])
    obs20, obst20 = tb.windowMean(Bdoth, time=smstat['time'], winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)


    obs_event = np.asarray(obs20) >= threshold

    minlen = min(len(run20), len(obs20))
    predicted_event = predicted_event[:minlen]
    obs_event = obs_event[:minlen]

    return predicted_event, obs_event


def write_table(ctable, date, stat, keywrd, threshold):
    f = open('ctables/{0}_{1}_ctable_{2}.txt'.format(date,stat, keywrd), 'a')
    f.write('\n')
    f.write('==============================\n')
    f.write('Threshold: {}\n'.format(threshold))
    f.write('==============================\n')
    f.write('\n')
    f.close()
    sys.stdout = open('ctables/{0}_{1}_ctable_{2}.txt'.format(date,stat, keywrd), 'a')
    ctable.summary(ci='bootstrap', verbose=True)
    sys.stdout = sys.__stdout__

def make_table(smdata, stations, outputdir, date, starttime, thresholds, keywrd):
    for stat in stations:
        dBdth_total, simtime_total = read_output_data(outputdir, keywrd, stat, starttime)
        smstat = smdata.station[stat]
        for threshold in thresholds:
            predicted_event, obs_event = make_boolean_arrays(dBdth_total, smstat, simtime_total, starttime, threshold)
            ctable = verify.Contingency2x2.fromBoolean(predicted_event, obs_event)

            print(len(np.where(predicted_event == True)[0]), threshold, stat)
            print(len(np.where(predicted_event == False)[0]), threshold, stat)
            print(simtime_total[0], starttime, simtime_total[-1])
            write_table(ctable, date, stat, keywrd, threshold)

if __name__ == "__main__":
    main(sys.argv)