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
import matplotlib.colors as colors
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
    #thresholds = [0.01]

    #make_table(smdata,stations,outputdir,date,starttime, thresholds,'unsmoothed')
    #make_table(smdata,stations,outputdir,date,starttime, thresholds,'hour')
    #make_table(smdata,stations,outputdir,date,starttime, thresholds,'30min')

    #midlat_sats = ['OTT', 'NEW', 'WNG', 'MEA']
    midlat_sats = ['BEL', 'CLF', 'FMC', 'HAD', 'MEA', 'OTT', 'SIT', 'THY', 'WNG', 'DOU', 'FUR', 'HLP', 'PIN', 'STJ', 'UPS', 'BFE', 'ESK', 'GIM', 'NEW', 'PBQ', 'SUA', 'VAL', 'FCC', 'IRT', 'NGK', 'RAL', 'TAR', 'VIC']

    highlat_sats= ['ABK', 'BLC', 'BRW', 'BJN', 'CBB', 'CMO', 'DNB', 'DOB', 'EAG','FSP','SMI','HRN','IQA','STF','KEV','KUV','LER','LYR','NAQ','NAL','NRD','NUR','OUJ','THL','RAN','RES','SVS','TAL','AMK','TIK','YKC']
    #highlat_sats = ['ABK', 'YKC', 'IQA']
    print('high lat')
    #grouping(outputdir, smdata, thresholds, highlat_sats, date, starttime)
    print('midlat')
    #grouping(outputdir, smdata, thresholds, midlat_sats, date, starttime)

    starttimes = ['20061214120000','20010831000000','20050831100000','20100405000000','20110805090000']

    stations = midlat_sats + highlat_sats
    #cross_event_grouping(outputdir, thresholds, highlat_sats, starttimes)
    #cross_event_grouping(outputdir, thresholds, midlat_sats, starttimes)
    #grouping(outputdir, smdata, thresholds, midlat_sats, 'Mid_Latitude', date, starttime)
    #grouping(outputdir, smdata, thresholds, highlat_sats, 'High_Latitude', date, starttime)
    create_polarplot(thresholds, stations, starttimes)


def grouping(outputdir, smdata, thresholds, stations, date, starttime):
    for threshold in thresholds:
        for keywrd in ['unsmoothed','30min','hour']:
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
                print(threshold, keywrd, starttime)
                date = starttime[:8]
                starttime = dt.datetime.strptime(starttime, '%Y%m%d%H%M%S')
                if date == '20010831': date = '20010830'
                outputdir = './outputs/{}/'.format(date)
                blockPrint()
                smdata = supermag_parser.supermag_parser('./supermag_data/{0}-supermag.txt'.format(date))
                enablePrint()
                origdata = open_output_data(outputdir, keywrd, starttime)
                for stat in stations:
                    if stat not in smdata.station: print(stat, 'not found')
                    else:
                        dBdth_total, simtime_total = just_read_output_data(origdata, keywrd, stat, starttime)
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
def create_dics(thresholds, stations, starttimes):
    dic = {}
    midlat_sats = ['BEL', 'CLF', 'FMC', 'HAD', 'MEA', 'OTT', 'SIT', 'THY', 'WNG', 'DOU', 'FUR', 'HLP', 'PIN', 'STJ', 'UPS', 'BFE', 'ESK', 'GIM', 'NEW', 'PBQ', 'SUA', 'VAL', 'FCC', 'IRT', 'NGK', 'RAL', 'TAR', 'VIC']
    highlat_sats= ['ABK', 'BLC', 'BRW', 'BJN', 'CBB', 'CMO', 'DNB', 'DOB', 'EAG','FSP','SMI','HRN','IQA','STF','KEV','KUV','LER','LYR','NAQ','NAL','NRD','NUR','OUJ','THL','RAN','RES','SVS','TAL','AMK','TIK','YKC']

    print('Initializing Dictionary')
    for threshold in thresholds:
        dic[threshold] = {}
        for keywrd in ['hour', 'unsmoothed', '30min']:
            dic[threshold][keywrd] = {}
            dic[threshold][keywrd]['highlat'] = {}
            dic[threshold][keywrd]['midlat'] = {}
            for i in range(5)[1:]:
                dic[threshold][keywrd]['highlat'][str(i)] = {}
                dic[threshold][keywrd]['midlat'][str(i)] = {}

                dic[threshold][keywrd]['highlat'][str(i)]['obs'] = []
                dic[threshold][keywrd]['midlat'][str(i)]['obs'] = []

                dic[threshold][keywrd]['highlat'][str(i)]['predicted'] = []
                dic[threshold][keywrd]['midlat'][str(i)]['predicted'] = []
    proc_dic={}
    for starttime in starttimes:
        proc_dic[starttime] = {}

        for station in stations:
            proc_dic[starttime][station] = {}
            for keywrd in ['hour', 'unsmoothed', '30min']:
                proc_dic[starttime][station][keywrd] = {}
                for threshold in thresholds:
                    proc_dic[starttime][station][keywrd][threshold] = {}  

    print('Filling Processing Dictionary')
    for keywrd in ['hour', 'unsmoothed', '30min']:
        for starttime in starttimes:
            strstart = starttime
            date = starttime[:8]
            print('Date: ', date)
            starttime = dt.datetime.strptime(starttime, '%Y%m%d%H%M%S')
            if date == '20010831': date = '20010830'
            outputdir = './outputs/{}/'.format(date)

            blockPrint()
            smdata = supermag_parser.supermag_parser('./supermag_data/{0}-supermag.txt'.format(date))
            enablePrint()
            origdata = open_output_data(outputdir, keywrd, starttime)
            for stat in stations:
                if stat not in smdata.station: continue
                else:                
                    smstat = smdata.station[stat]
                    dBdth, simtime = process_dBdth(origdata, stat)
                    mlt_len = len(smstat['mlt'])
                    stat_mlt20 = tb.windowMean(smstat['mlt'], time=smstat['time'][:mlt_len], winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
                    for threshold in thresholds:
                        predicted_event_sat, obs_event_sat = make_boolean_arrays(dBdth, smstat, simtime, starttime, threshold)   
                        proc_dic[strstart][station][keywrd][threshold]['predicted_events'] = predicted_event_sat
                        proc_dic[strstart][station][keywrd][threshold]['obs_events'] = obs_event_sat

                        proc_dic[strstart][station][keywrd][threshold]['mlt'] = stat_mlt20[0]
                        proc_dic[strstart][station][keywrd][threshold]['mlat'] = smstat['mlat']

    # proc_dic is filled with boolean arrays for each threshold level, event, and satellite
    # next, we need to seciton up the data into mlt sectors
    

    print('Filling Main Dictionary')

    for threshold in thresholds:
        for keywrd in ['hour', 'unsmoothed', '30min']:
            for stat in stations:
                if stat in highlat_sats or stat in midlat_sats:
                    #print(np.where(proc_dic[strstart][station][keywrd][threshold]['obs_events'] == True))
                    #print(np.where(proc_dic[strstart][station][keywrd][threshold]['predicted_events'] == True))
                    for i in range(len(proc_dic[strstart][station][keywrd][threshold]['mlt'])):
                        curmlt = proc_dic[strstart][station][keywrd][threshold]['mlt'][i]
                        print(curmlt)
                        flag = ''
                        if stat in highlat_sats:
                            flag = 'highlat'
                        else:
                            flag = 'midlat'


                        if 0 <= curmlt < 6:
                            dic[threshold][keywrd][flag]['1']['obs'] += [proc_dic[strstart][station][keywrd][threshold]['obs_events'][i]]
                            dic[threshold][keywrd][flag]['1']['predicted'] += [proc_dic[strstart][station][keywrd][threshold]['predicted_events'][i]]

                        elif 6 <= curmlt < 12:
                            
                            dic[threshold][keywrd][flag]['2']['obs'] += [proc_dic[strstart][station][keywrd][threshold]['obs_events'][i]]
                            dic[threshold][keywrd][flag]['2']['predicted'] += [proc_dic[strstart][station][keywrd][threshold]['predicted_events'][i]]                            

                        elif 12<= curmlt < 18:
                            
                            dic[threshold][keywrd][flag]['3']['obs'] += [proc_dic[strstart][station][keywrd][threshold]['obs_events'][i]]
                            dic[threshold][keywrd][flag]['3']['predicted'] += [proc_dic[strstart][station][keywrd][threshold]['predicted_events'][i]]

                        elif 18<= curmlt < 24:

                            dic[threshold][keywrd][flag]['4']['obs'] += [proc_dic[strstart][station][keywrd][threshold]['obs_events'][i]]
                            dic[threshold][keywrd][flag]['4']['predicted'] += [proc_dic[strstart][station][keywrd][threshold]['predicted_events'][i]]

    
    # now calculate tt for each sector to get metrics
    print('lets check differences')
    print(dic[0.3]['unsmoothed']['highlat']['1']['obs'] == dic[0.3]['unsmoothed']['midlat']['1']['obs'])
    print(dic[0.3]['unsmoothed']['highlat']['2']['obs'] == dic[0.3]['unsmoothed']['midlat']['2']['obs'])
    print(dic[0.3]['unsmoothed']['highlat']['3']['obs'] == dic[0.3]['unsmoothed']['midlat']['3']['obs'])
    print('calculating statistics')
    for threshold in thresholds:
        for keywrd in ['hour', 'unsmoothed', '30min']:
            for i in range(5)[1:]:
                for flag in ['highlat', 'midlat']:
                    obs_events = dic[threshold][keywrd][flag][str(i)]['obs']
                    predic_events = dic[threshold][keywrd][flag][str(i)]['predicted']
                    tt = verify.Contingency2x2.fromBoolean(predic_events, obs_events)

                    sector_skill = tt.heidke(ci='bootstrap')
                    dic[threshold][keywrd][flag][str(i)]['heidke'] = sector_skill[0]

 

    return proc_dic, dic

def create_polarplot(thresholds, stations, starttimes):
    
    proc_dic, dic = create_dics(thresholds, stations, starttimes)
    
    midlat_sats = ['BEL', 'CLF', 'FMC', 'HAD', 'MEA', 'OTT', 'SIT', 'THY', 'WNG', 'DOU', 'FUR', 'HLP', 'PIN', 'STJ', 'UPS', 'BFE', 'ESK', 'GIM', 'NEW', 'PBQ', 'SUA', 'VAL', 'FCC', 'IRT', 'NGK', 'RAL', 'TAR', 'VIC']
    highlat_sats= ['ABK', 'BLC', 'BRW', 'BJN', 'CBB', 'CMO', 'DNB', 'DOB', 'EAG','FSP','SMI','HRN','IQA','STF','KEV','KUV','LER','LYR','NAQ','NAL','NRD','NUR','OUJ','THL','RAN','RES','SVS','TAL','AMK','TIK','YKC']
    
    fig = plt.figure()
    ax = fig.add_subplot(111, polar = True)
    
    deltaskill_highlat = []
    for i in range(5)[1:]:
         delta = dic[0.3]['hour']['highlat'][str(i)]['heidke'] - dic[0.3]['unsmoothed']['highlat'][str(i)]['heidke']
         deltaskill_highlat += [delta]
    deltaskill_midlat = []
    for i in range(5)[1:]:
         delta =  dic[0.3]['hour']['midlat'][str(i)]['heidke'] - dic[0.3]['unsmoothed']['midlat'][str(i)]['heidke']
         deltaskill_midlat += [delta]
    
    
    #deltaskill_highlat = [np.NaN, -0.01591511936339518, -0.020682148040638018, -0.25263157894736854] 
    #deltaskill_midlat = [np.NaN, -0.01591511936339518, -0.020682148040638018, -0.25263157894736854]


    print(deltaskill_highlat, deltaskill_midlat)

    cdictpos = {'green':  ((0.0, 0.0, 0.0),   # no green at 0
                  (0.5, 0.4, 0.4),   # all channels set to 1.0 at 0.5 to create white
                  (1.0, 0.8, 0.8)),  # set to 0.8 so its not too bright at 1

        'red': ((0.0, 0.0, 0.0),   # set to 0.8 so its not too bright at 0
                  (0.5, 0.0, 0.0),   # all channels set to 1.0 at 0.5 to create white
                  (1.0, 0.0, 0.0)),  # no red at 1

        'blue':  ((1.0, 0.0, 0.0),   # no blue at 0
                  (0.5, 0.0, 0.0),   # all channels set to 1.0 at 0.5 to create white
                  (1.0, 0.0, 0.0))   # no blue at 1
       }
    cdictneg = {'green':  ((0.0, 0.0, 0.0),   # no green at 0
                  (0.5, 0.0, 0.0),   # all channels set to 1.0 at 0.5 to create white
                  (1.0, 0.0, 0.0)),  # set to 0.8 so its not too bright at 1

        'red': ((0.0, 0.0, 0.0),   # set to 0.8 so its not too bright at 0
                  (0.5, 0.4, 0.4),   # all channels set to 1.0 at 0.5 to create white
                  (1.0, 0.8, 0.8)),  # no red at 1

        'blue':  ((0.0, 0.0, 0.0),   # no blue at 0
                  (0.5, 0.0, 0.0),   # all channels set to 1.0 at 0.5 to create white
                  (1.0, 0.0, 0.0))   # no blue at 1
       }
    #GnRdPos = colors.LinearSegmentedColormap('GnRdPos', cdictpos)
    #GnRdNeg = colors.LinearSegmentedColormap('GnRdNeg', cdictneg)

    GnRdPos = plt.get_cmap('Greens')
    GnRdNeg = plt.get_cmap('Reds')

    norm = colors.Normalize(vmin=0, vmax=0.5)

    for i in range(4):
        if np.isnan(deltaskill_highlat[i]): color = 'white'
        elif deltaskill_highlat[i] < 0:
            color = GnRdNeg(norm(abs(deltaskill_highlat[i])))
        else: 
            color = GnRdPos(norm(deltaskill_highlat[i]))
        ax.bar(7*np.pi/4 + i*np.pi/2, 1, width=1 * np.pi / 2, bottom=0,
               color=color, edgecolor = color)
        
        # mid lat
        if np.isnan(deltaskill_midlat[i]): color = 'white'
        elif deltaskill_midlat[i] < 0:
            color = GnRdNeg(norm(abs(deltaskill_midlat[i])))
        else:
            color = GnRdPos(norm(deltaskill_midlat[i]))

        ax.bar(7*np.pi/4 + i*np.pi/2, 0.5, width=np.pi / 2, bottom=0,
               color=color, edgecolor = color)

    '''
    if np.isnan(deltaskill_highlat[1]): color = 'white'
    elif deltaskill_highlat[1] < 0:
        color = GnRdNeg(norm(abs(deltaskill_highlat[1])))
    else:
        color = GnRdPos(norm(deltaskill_highlat[1]))
    ax.bar(np.pi/4, 1, width=1 * np.pi / 2, bottom=0,
           color=color, edgecolor = color)

    if np.isnan(deltaskill_highlat[2]): color = 'white'
    elif deltaskill_highlat[2] < 0:
        color = GnRdNeg(norm(abs(deltaskill_highlat[2])))
    else:
        color = GnRdPos(norm(deltaskill_highlat[2]))

    ax.bar(3*np.pi/4, 1, width=1 * np.pi / 2, bottom=0,
           color=color, edgecolor = color)

    if np.isnan(deltaskill_highlat[3]): color = 'white'
    elif deltaskill_highlat[3] < 0:
        color = GnRdNeg(norm(abs(deltaskill_highlat[3])))
    else:
        color = GnRdPos(norm(deltaskill_highlat[3]))

    ax.bar(5*np.pi/4, 1, width=1 * np.pi / 2, bottom=0,
           color=color, edgecolor = color)

    
    if np.isnan(deltaskill_midlat[0]): color = 'white'
    elif deltaskill_midlat[0] < 0:
        color = GnRdNeg(norm(abs(deltaskill_midlat[0])))
    else:
        color = GnRdPos(norm(deltaskill_midlat[0]))

    ax.bar(7*np.pi/4, 0.5, width=np.pi / 2, bottom=0,
           color=color, edgecolor = color)


    if np.isnan(deltaskill_midlat[1]): color = 'white'
    elif deltaskill_midlat[1] < 0:
        color = GnRdNeg(norm(abs(deltaskill_midlat[1])))
    else:
        color = GnRdPos(norm(deltaskill_midlat[1]))

    ax.bar(np.pi/4, 0.5, width=np.pi / 2, bottom=0,
           color=color, edgecolor = color)


    if np.isnan(deltaskill_midlat[2]): color = 'white'
    elif deltaskill_midlat[2] < 0:
        color = GnRdNeg(norm(abs(deltaskill_midlat[2])))
    else:
        color = GnRdPos(norm(deltaskill_midlat[2]))

    ax.bar(3*np.pi/4, 0.5, width=np.pi / 2, bottom=0,
           color=color, edgecolor = color)

    if np.isnan(deltaskill_midlat[3]): color = 'white'
    elif deltaskill_midlat[3] < 0:
        color = GnRdNeg(norm(abs(deltaskill_midlat[3])))
    else:
        color = GnRdPos(norm(deltaskill_midlat[3]))

    ax.bar(5*np.pi/4, 0.5, width=np.pi / 2, bottom=0,
           color=color, edgecolor = color)
    '''

    plt.show()



                
    
def open_output_data(outputdir, keywrd, starttime):
    files = sorted(glob.glob(outputdir + '/{}/mag*mag'.format(keywrd)))
    for fpath in files:
        origdata = spacepy.pybats.bats.MagFile(fpath)
        origdata.calc_h()
        origdata.calc_dbdt() 
    return origdata

def process_dBdth(origdata, stat):
    dBdth_total = []
    simtime_total = []
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

def read_output_data(outputdir, keywrd, stat, starttime):
    origdata = open_output_data(outputdir, keywrd, starttime)
    dBdth_total, simtime_total = process_dBdth(origdata, stat)

    return dBdth_total, simtime_total

def just_read_output_data(origdata, keywrd, stat, starttime):
    
    dBdth_total, simtime_total = process_dBdth(origdata, stat)

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