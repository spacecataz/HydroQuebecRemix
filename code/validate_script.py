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
import matplotlib

# Example command to run: python validate_script.py 20061214120000
# As it is right now, its only doing the polar plot so this command will run it!

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

def main(args):
    # there are two ways to run this, one specifying both the date and the start time, or just the start time
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

    #midlat_stats = ['OTT', 'NEW', 'WNG', 'MEA']
    midlat_stats = ['BEL', 'CLF', 'FMC', 'HAD', 'MEA', 'OTT', 'SIT', 'THY', 'WNG', 'DOU', 'FUR', 'HLP', 'PIN', 'STJ', 'UPS', 'BFE', 'ESK', 'GIM', 'NEW', 'PBQ', 'SUA', 'VAL', 'FCC', 'IRT', 'NGK', 'RAL', 'TAR', 'VIC']

    highlat_stats= ['ABK', 'BLC', 'BRW', 'BJN', 'CBB', 'CMO', 'DNB', 'DOB', 'EAG','FSP','SMI','HRN','IQA','STF','KEV','KUV','LER','LYR','NAQ','NAL','NRD','NUR','OUJ','THL','RAN','RES','SVS','TAL','AMK','TIK','YKC']
    #highlat_stats = ['ABK', 'YKC', 'IQA']

    # this seciton does grouping by high lat and mid lat for single event
    print('high lat')
    #grouping(outputdir, smdata, thresholds, highlat_stats, date, starttime)
    print('midlat')
    #grouping(outputdir, smdata, thresholds, midlat_stats, date, starttime)

    # this section does grouping across all events for mid lat and high lat
    starttimes = ['20061214120000','20010831000000','20050831100000','20100405000000','20110805090000']
    #starttimes = ['20061214120000','20010831000000']
    stations = midlat_stats + highlat_stats
    #cross_event_grouping(outputdir, thresholds, highlat_stats, starttimes)
    #cross_event_grouping(outputdir, thresholds, midlat_stats, starttimes)
    #grouping(outputdir, smdata, thresholds, midlat_stats, 'Mid_Latitude', date, starttime)
    #grouping(outputdir, smdata, thresholds, highlat_stats, 'High_Latitude', date, starttime)
    
    # this creates the polar plot (or tries to...)
    create_polarplot(thresholds, stations, starttimes)


def grouping(outputdir, smdata, thresholds, stations, date, starttime):
    # This function groups across stations, specified in the input 'stations'
    for threshold in thresholds:
        for keywrd in ['unsmoothed','30min','hour']:
            predicted_event_tot = []
            obs_event_tot = []
            for stat in stations:
                if stat not in smdata.station: continue
                else:
                    dBdth_total, simtime_total = read_output_data(outputdir, keywrd, stat, starttime)
                    smstat = smdata.station[stat]
                    predicted_event_stat, obs_event_stat = make_boolean_arrays(dBdth_total, smstat, simtime_total, starttime, threshold)
                    predicted_event_tot += predicted_event_stat.tolist()
                    obs_event_tot += obs_event_stat.tolist()

            ctable = verify.Contingency2x2.fromBoolean(predicted_event_tot, obs_event_tot)
            if stations[0] == 'ABK':
                group = 'highlat'
            else:
                group = 'midlat'
            write_table(ctable, date, group, keywrd, threshold)

def cross_event_grouping(outputdir, thresholds, stations, starttimes):
    # This function will group across events for a given set of statations
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
                        predicted_event_stat, obs_event_stat = make_boolean_arrays(dBdth_total, smstat, simtime_total, starttime, threshold)
                        predicted_event_tot += predicted_event_stat.tolist()
                        obs_event_tot += obs_event_stat.tolist()
            print('Lengths before CTable')
            print(len(predicted_event_tot), len(obs_event_tot))
            ctable = verify.Contingency2x2.fromBoolean(predicted_event_tot, obs_event_tot)
            if stations[0] == 'ABK':
                group = 'combined_highlat'
            else:
                group = 'combined_midlat'
            write_table(ctable, date, group, keywrd, threshold)

def create_dics(thresholds, stations, starttimes):
    # This is a helper funciton for polar plotting
    # function creates the dictionaries used to perform the statistics

    dic = {}
    midlat_stats = ['BEL', 'CLF', 'FMC', 'HAD', 'MEA', 'OTT', 'SIT', 'THY', 'WNG', 'DOU', 'FUR', 'HLP', 'PIN', 'STJ', 'UPS', 'BFE', 'ESK', 'GIM', 'NEW', 'PBQ', 'SUA', 'VAL', 'FCC', 'IRT', 'NGK', 'RAL', 'TAR', 'VIC']
    highlat_stats= ['ABK', 'BLC', 'BRW', 'BJN', 'CBB', 'CMO', 'DNB', 'DOB', 'EAG','FSP','SMI','HRN','IQA','STF','KEV','KUV','LER','LYR','NAQ','NAL','NRD','NUR','OUJ','THL','RAN','RES','SVS','TAL','AMK','TIK','YKC']

    print('Initializing Dictionary')
    # This dicitonary is organized by Threshold Value/Smoothing Level/High or Mid Lat/Sector/Predicted or Observed
    # Stores boolean values
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
    # This dictionary is used to first load all of the data
    # It stores the boolean values 
    # Organized by Date/Station/Smoothing Level/Threshold
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
    # This dictionary is organized by date/station/smoothing level/threshold
    # it is filled with boolean values (based off of each threshold)
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
                if stat not in smdata.station: 
                    if stat in proc_dic[strstart].keys(): del proc_dic[strstart][stat]
                else: 
                    # Get boolean values           
                    smstat = smdata.station[stat]
                    dBdth, simtime = process_dBdth(origdata, stat)
                    mlt_len = len(smstat['mlt'])
                    stat_mlt20 = tb.windowMean(smstat['mlt'], time=smstat['time'][:mlt_len], winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime)
                    for threshold in thresholds:
                        predicted_event_stat, obs_event_stat = make_boolean_arrays(dBdth, smstat, simtime, starttime, threshold)   

                        proc_dic[strstart][stat][keywrd][threshold]['predicted'] = predicted_event_stat
                        proc_dic[strstart][stat][keywrd][threshold]['obs'] = obs_event_stat

                        proc_dic[strstart][stat][keywrd][threshold]['mlt'] = stat_mlt20[0]
                        proc_dic[strstart][stat][keywrd][threshold]['mlat'] = smstat['mlat']

    # proc_dic is filled with boolean arrays for each threshold level, event, and station
    # next, we need to seciton up the data into mlt sectors
    print('Filling Main Dictionary')
    # This dictionary is organized 
    for threshold in thresholds:
        for keywrd in ['hour', 'unsmoothed', '30min']:
            for stat in stations:
                if stat in highlat_stats or stat in midlat_stats:
                    flag = ''
                    if stat in highlat_stats:
                        flag = 'highlat'
                    else:
                        flag = 'midlat'  
                    if stat not in proc_dic[strstart].keys(): continue
                    else:
                        # '''
                        mlt = np.array(proc_dic[strstart][stat][keywrd][threshold]['mlt'])
                        sector1 = np.logical_or(mlt >= 21, mlt < 3)
                        sector2 = np.logical_and(mlt >= 3, mlt < 9)
                        sector3 = np.logical_and(mlt >= 9, mlt <15)
                        sector4 = np.logical_and(mlt >= 15, mlt < 21)

                        for k in ['obs', 'predicted']:
                            dic[threshold][keywrd][flag]['1'][k] += list(np.array(proc_dic[strstart][stat][keywrd][threshold][k])[sector1])
                            dic[threshold][keywrd][flag]['2'][k] += list(np.array(proc_dic[strstart][stat][keywrd][threshold][k])[sector2])
                            dic[threshold][keywrd][flag]['3'][k] += list(np.array(proc_dic[strstart][stat][keywrd][threshold][k])[sector3])
                            dic[threshold][keywrd][flag]['4'][k] += list(np.array(proc_dic[strstart][stat][keywrd][threshold][k])[sector4])

    # now calculate tt for each sector to get metrics
    print('calculating statistics')
    for threshold in thresholds:
        for keywrd in ['hour', 'unsmoothed', '30min']:
            for i in range(1,5):
                for flag in ['highlat', 'midlat']:
                    obs_events = dic[threshold][keywrd][flag][str(i)]['obs']
                    predic_events = dic[threshold][keywrd][flag][str(i)]['predicted']

                    # These are useful for tallying events
                    #print(threshold, keywrd, i, flag, sum(obs_events))
                    #print(sum(predic_events))
                    tt = verify.Contingency2x2.fromBoolean(predic_events, obs_events)

                    sector_skill = tt.heidke(ci='bootstrap')
                    dic[threshold][keywrd][flag][str(i)]['heidke'] = sector_skill[0]

 

    return proc_dic, dic

def create_polarplot(thresholds, stations, starttimes):

    proc_dic, dic = create_dics(thresholds, stations, starttimes)
    
    midlat_stats = ['BEL', 'CLF', 'FMC', 'HAD', 'MEA', 'OTT', 'SIT', 'THY', 'WNG', 'DOU', 'FUR', 'HLP', 'PIN', 'STJ', 'UPS', 'BFE', 'ESK', 'GIM', 'NEW', 'PBQ', 'SUA', 'VAL', 'FCC', 'IRT', 'NGK', 'RAL', 'TAR', 'VIC']
    highlat_stats= ['ABK', 'BLC', 'BRW', 'BJN', 'CBB', 'CMO', 'DNB', 'DOB', 'EAG','FSP','SMI','HRN','IQA','STF','KEV','KUV','LER','LYR','NAQ','NAL','NRD','NUR','OUJ','THL','RAN','RES','SVS','TAL','AMK','TIK','YKC']


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
    

    print(deltaskill_highlat, deltaskill_midlat)

    #GnRdPos = plt.get_cmap('Greens')
    #GnRdNeg = plt.get_cmap('Reds')
    cdict = {'red':  ((0.0, 0.8, 0.8),   # red at 0
                      (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                      (1.0, 0.0, 0.0)),  # no red at 1

            'green': ((0.0, 0.0, 0.0),   # no green at 0
                      (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                      (1.0, 0.8, 0.8)),  # green at 1

            'blue':  ((0.0, 0.0, 0.0),   # no blue at 0
                      (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                      (1.0, 0.0, 0.0))   # no blue at 1
           }
    # Create the colormap
    GnRd = colors.LinearSegmentedColormap('GnRd', cdict)
    # Create the normalization (without norm only 0-1 values are allowed)
    norm = colors.Normalize(vmin=-0.5, vmax=0.5)

    # Plot each sector
    for i in range(4):
        # Plot high lat
        color = GnRd(norm(deltaskill_highlat[i]))
        ax.bar(i*np.pi/2, 1, width=1 * np.pi / 2, bottom=0,
               color=color, edgecolor = color)
        
        # Plot mid lat
        color = GnRd(norm(deltaskill_midlat[i]))
        ax.bar(i*np.pi/2, 0.5, width=np.pi / 2, bottom=0,
               color=color, edgecolor = color)

    # Add in colorbar
    clb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=GnRd), ax=ax)
    clb.ax.set_title(r'$\Delta$ Skill')
    ax.set_thetagrids((270, 0, 90, 180), labels = ['00','6', '12', '18'])
    ax.set_yticklabels([])
    fig.suptitle('Comparing Change in Skill Between Hourly and Unsmoothed Inputs')
    fig.savefig('hourly_polarplot.png')

    

    fig = plt.figure()
    ax = fig.add_subplot(111, polar = True)

    deltaskill_highlat = []
    for i in range(5)[1:]:
         delta = dic[0.3]['30min']['highlat'][str(i)]['heidke'] - dic[0.3]['unsmoothed']['highlat'][str(i)]['heidke']
         deltaskill_highlat += [delta]
    deltaskill_midlat = []
    for i in range(5)[1:]:
         delta =  dic[0.3]['30min']['midlat'][str(i)]['heidke'] - dic[0.3]['unsmoothed']['midlat'][str(i)]['heidke']
         deltaskill_midlat += [delta]

    #print(deltaskill_highlat, deltaskill_midlat)


    for i in range(4):
        color = GnRd(norm(deltaskill_highlat[i]))
        ax.bar(i*np.pi/2, 1, width=1 * np.pi / 2, bottom=0,
               color=color, edgecolor = color)
        
        # mid lat
        color = GnRd(norm(deltaskill_midlat[i]))
        ax.bar(i*np.pi/2, 0.5, width=np.pi / 2, bottom=0,
               color=color, edgecolor = color)

    clb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=GnRd), ax=ax)
    clb.ax.set_title(r'$\Delta$ Skill')
    ax.set_thetagrids((270, 0, 90, 180), labels = ['00','6', '12', '18'])
    fig.suptitle('Comparing Change in Skill Between 30min and Unsmoothed Inputs')
    ax.set_yticklabels([])
    fig.savefig('30min_polarplot.png')
    plt.show()
    

    return fig, ax
    
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
    #dBdth = np.array([linalg.norm([dBdtn[i], dBdte[i]], axis = 0)])
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

            write_table(ctable, date, stat, keywrd, threshold)

if __name__ == "__main__":
    midlat_stats = ['BEL', 'CLF', 'FMC', 'HAD', 'MEA', 'OTT', 'SIT', 'THY', 'WNG', 'DOU', 'FUR', 'HLP', 'PIN', 'STJ', 'UPS', 'BFE', 'ESK', 'GIM', 'NEW', 'PBQ', 'SUA', 'VAL', 'FCC', 'IRT', 'NGK', 'RAL', 'TAR', 'VIC']

    highlat_stats= ['ABK', 'BLC', 'BRW', 'BJN', 'CBB', 'CMO', 'DNB', 'DOB', 'EAG','FSP','SMI','HRN','IQA','STF','KEV','KUV','LER','LYR','NAQ','NAL','NRD','NUR','OUJ','THL','RAN','RES','SVS','TAL','AMK','TIK','YKC']
    #highlat_stats = ['ABK', 'YKC', 'IQA']

    # this seciton does grouping by high lat and mid lat for single event
    print('high lat')
    #grouping(outputdir, smdata, thresholds, highlat_stats, date, starttime)
    print('midlat')
    #grouping(outputdir, smdata, thresholds, midlat_stats, date, starttime)

    # this section does grouping across all events for mid lat and high lat
    starttimes = ['20061214120000','20010831000000','20050831100000','20100405000000','20110805090000']
    #starttimes = ['20061214120000','20010831000000']
    stations = midlat_stats + highlat_stats
    thresholds = [0.3, 0.7, 1.1, 1.5] #nT/s
    fig, ax = create_polarplot(thresholds, stations, starttimes)

    #main(sys.argv)