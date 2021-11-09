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

    #make_table(smdata,stations,outputdir,date,starttime, thresholds,'1min')
    #make_table(smdata,stations,outputdir,date,starttime, thresholds,'60min')
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
        for keywrd in ['1min','15min','30min','60min']:
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

def cross_event_grouping(thresholds, stations, starttimes):
    # This function will group across events for a given set of statations
    for threshold in thresholds:
        for keywrd in ['1min','15min','30min','60min']:   
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
        for keywrd in ['1min','15min','30min','60min']:
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
            for keywrd in ['1min','15min','30min','60min','obs']:
                proc_dic[starttime][station][keywrd] = {}
                for threshold in thresholds:
                    proc_dic[starttime][station][keywrd][threshold] = {}  

    print('Filling Processing Dictionary')
    # This dictionary is organized by date/station/smoothing level/threshold
    # it is filled with boolean values (based off of each threshold)
    for keywrd in ['1min','15min','30min','60min']:
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
                    if stat in proc_dic[strstart].keys(): 
                        del proc_dic[strstart][stat]
                        print(stat, 'not in SM Data', starttime)
                else: 
                    # Get boolean values           
                    smstat = smdata.station[stat]
                    dBdth, simtime = process_dBdth(origdata, stat)

                    mlt_len = len(smstat['mlt'])
                    stat_mlt20 = tb.windowMean(smstat['mlt'], time=smstat['time'][:mlt_len], winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime)


                    for threshold in thresholds:
                        predicted_event_stat, obs_event_stat = make_boolean_arrays(dBdth, smstat, simtime, starttime, threshold)   

                        proc_dic[strstart][stat][keywrd][threshold]['predicted'] = predicted_event_stat
                        if threshold == 0.3:

                            subset = origdata[stat]
                            dBe = decimate(subset['dBe'], 6)
                            dBn = decimate(subset['dBn'], 6)
                            Bh = np.array([linalg.norm([dBn[i], dBe[i]]) for i in range(len(dBn))])

                            mltmin, mlttmin = tb.windowMean(smstat['mlt'], time=smstat['time'], winsize=dt.timedelta(minutes=1), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
                            datamin, datatmin = tb.windowMean(Bh, time=subset['time'][::6], winsize=dt.timedelta(minutes=1), overlap=dt.timedelta(0), st_time=starttime, op=np.max)

                            Bdoth = np.array([linalg.norm(smstat['B'][i,:2]) for i in range(len(smstat['B']))])
                            obs20, obst20 = tb.windowMean(Bdoth, time=smstat['time'], winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
                            
                            minlen = min(len(mltmin), len(datamin))
                            mltmin = mltmin[:minlen]
                            datamin = datamin[:minlen]
                            obs20 = obs20[:minlen]

                            proc_dic[strstart][stat][keywrd][threshold]['predicted_data'] = datamin
                            proc_dic[strstart][stat]['obs'][threshold]['predicted_data'] = obs20
                            proc_dic[strstart][stat]['obs'][threshold]['predicted_mlt'] = mltmin
                            proc_dic[strstart][stat][keywrd][threshold]['predicted_mlt'] = mltmin

                        proc_dic[strstart][stat][keywrd][threshold]['obs'] = obs_event_stat

                        proc_dic[strstart][stat][keywrd][threshold]['mlt'] = stat_mlt20[0]
                        proc_dic[strstart][stat][keywrd][threshold]['mlat'] = smstat['mlat']

    # proc_dic is filled with boolean arrays for each threshold level, event, and station
    # next, we need to seciton up the data into mlt sectors
    print('Filling Main Dictionary')
    # This dictionary is organized 
    for threshold in thresholds:
        for keywrd in ['1min','15min','30min','60min']:
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
                 
    # power spectra
    print('calculating power spectra')

    for starttime in starttimes:
        power_spectra_plots(starttime,stations,proc_dic,'all')
        power_spectra_plots(starttime,highlat_stats,proc_dic,'high_lat')
        power_spectra_plots(starttime,midlat_stats,proc_dic,'mid_lat')

    return

    # now calculate tt for each sector to get metrics
    print('calculating statistics')
    for threshold in thresholds:
        for keywrd in ['1min','15min','30min','60min']:
            for i in range(1,5):
                for flag in ['highlat', 'midlat']:
                    obs_events = dic[threshold][keywrd][flag][str(i)]['obs']
                    predic_events = dic[threshold][keywrd][flag][str(i)]['predicted']

                    # These are useful for tallying events
                    #print(threshold, keywrd, i, flag, sum(obs_events))
                    #print(sum(predic_events))
                    tt = verify.Contingency2x2.fromBoolean(predic_events, obs_events)


                    sector_skill = tt.heidke(ci='bootstrap')
                    sector_bias = tt.bias(ci='bootstrap')
                    write_table(tt, 'sector', i, keywrd, threshold)
                    dic[threshold][keywrd][flag][str(i)]['heidke'] = sector_skill[0]
                    dic[threshold][keywrd][flag][str(i)]['bias'] = sector_bias[0]

 

    return proc_dic, dic

def sector_pad(data,sector_array):
    if len(sector_array) > len(data):
        sectori_data = data
    else:
        sectori_data = np.array(data)[sector_array]

    diff = 1000-len(sectori_data)
    diffmod = diff%2
    begin = np.zeros(int(diff/2))
    end = np.zeros(int(diff/2) + diffmod)
    sectori_data = np.append(begin,sectori_data)
    sectori_data = np.append(sectori_data,end)              
    return sectori_data

def power_spectra_plots(starttime,stations,proc_dic,grouping):
    total_power = []
    total_power_1 = []
    total_power_2 = []
    total_power_3 = []
    total_power_4 = []

    power_count = 0

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    fig4, ax4 = plt.subplots()
    fig5, ax5 = plt.subplots()

    for stat in stations:
        if stat in proc_dic[starttime].keys(): 
            sector1_list= []
            sector2_list=[]
            sector3_list=[]
            sector4_list=[]
            data_list=[]

            for keywrd in ['obs','1min','15min','30min','60min']:
                date = starttime[:8]
                mlt = np.array(proc_dic[starttime][stat][keywrd][0.3]['predicted_mlt'])
                data = np.array(proc_dic[starttime][stat][keywrd][0.3]['predicted_data'])
                data_list += [(keywrd,data)]

                sector1 = np.logical_or(mlt >= 21, mlt < 3)
                sector2 = np.logical_and(mlt >= 3, mlt < 9)
                sector3 = np.logical_and(mlt >= 9, mlt <15)
                sector4 = np.logical_and(mlt >= 15, mlt < 21)

                sector1_data = sector_pad(data,sector1)             
                sector1_list += [(keywrd,list(sector1_data))]

                sector2_data = sector_pad(data,sector2)             
                sector2_list += [(keywrd,list(sector2_data))]

                sector3_data = sector_pad(data,sector3)             
                sector3_list += [(keywrd,list(sector3_data))]

                sector4_data = sector_pad(data,sector4)             
                sector4_list += [(keywrd,list(sector4_data))]

            freq1,idx1 = power_spectra(sector1_list, stat, '1','combined',date)
            freq2, idx2 = power_spectra(sector2_list, stat, '2','combined',date)
            freq3, idx3 = power_spectra(sector3_list, stat, '3','combined',date)
            freq4, idx4 = power_spectra(sector4_list, stat, '4','combined',date)

            freq, idx = power_spectra(data_list, stat, 'total_data', 'combined', date)

            smoothing_levels= ['obs','1min','15min','30min','60min']
            for i in range(1,5):
                ax1.loglog(freq[i], idx[i],color='#a9a9a9')
                ax2.loglog(freq1[i], idx1[i],color='#a9a9a9')
                ax3.loglog(freq2[i], idx2[i],color='#a9a9a9')
                ax4.loglog(freq3[i], idx3[i],color='#a9a9a9')
                ax5.loglog(freq4[i], idx4[i],color='#a9a9a9')

            if total_power == []:
                total_power = idx
                total_power_1 = idx1
                total_power_2 = idx2
                total_power_3 = idx3
                total_power_4 = idx4
                power_count += 1

            else:
                for i in range(5):

                    total_power[i] += idx[i]

                    total_power_1[i] += idx1[i]
                    total_power_2[i] += idx2[i]
                    total_power_3[i] += idx3[i]
                    total_power_4[i] += idx4[i]
                    power_count += 1

    for i in range(5):
        total_power[i] = total_power[i]/power_count
        total_power_1[i] = total_power_1[i]/power_count
        total_power_2[i] = total_power_2[i]/power_count
        total_power_3[i] = total_power_3[i]/power_count
        total_power_4[i] = total_power_4[i]/power_count

    for i in range(1,5):
        ax1.loglog(freq[i], total_power[i], label=smoothing_levels[i])
        ax2.loglog(freq1[i], total_power_1[i], label=smoothing_levels[i])
        ax3.loglog(freq2[i], total_power_2[i], label=smoothing_levels[i])
        ax4.loglog(freq3[i], total_power_3[i], label=smoothing_levels[i])
        ax5.loglog(freq4[i], total_power_4[i], label=smoothing_levels[i])

    ax1.loglog(freq[0], total_power[0], 'k', label='Observed')
    ax2.loglog(freq1[0], total_power_1[0],'k', label='Observed')
    ax3.loglog(freq2[0], total_power_2[0],'k', label='Observed')
    ax4.loglog(freq3[0], total_power_3[0],'k', label='Observed')
    ax5.loglog(freq4[0], total_power_4[0],'k', label='Observed')

    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax5.legend()

    plt.ylabel(r'$[nT s]^2$')
    ax1.set_xlabel('[Hz]')
    ax2.set_xlabel('[Hz]')
    ax3.set_xlabel('[Hz]')
    ax4.set_xlabel('[Hz]')
    ax5.set_xlabel('[Hz]')

    ax1.set_title(date + grouping)
    ax2.set_title(date + ' Sector 1'+grouping)
    ax3.set_title(date+ ' Sector 2'+grouping)
    ax4.set_title(date+ ' Sector 3'+grouping)
    ax5.set_title(date+ ' Sector 4'+grouping)

    for ax in [ax1,ax2,ax3,ax4,ax5]:
        ax.axvline(x=0.001111,linestyle='--', color='black')

        ax.annotate("15min", xy=[0.001111, 0.4], fontsize=10,rotation=90)
        ax.annotate("30min", xy=[0.0005556, 0.4], fontsize=10,rotation=90)
        ax.annotate("60min", xy=[0.0002778, 0.4], fontsize=10,rotation=90)


        ax.axvline(x=0.0005556,linestyle='--', color='black')
        ax.axvline(x=0.0002778,linestyle='--', color='black')

    fig1.savefig('plots/powerspectra/{0}_combined_{1}_average_powerspectra.png'.format(date,grouping))
    fig2.savefig('plots/powerspectra/{0}_sector1_{1}_average_powerspectra.png'.format(date,grouping))
    fig3.savefig('plots/powerspectra/{0}_sector2_{1}_average_powerspectra.png'.format(date,grouping))
    fig4.savefig('plots/powerspectra/{0}_sector3_{1}_average_powerspectra.png'.format(date,grouping))
    fig5.savefig('plots/powerspectra/{0}_sector4_{1}_average_powerspectra.png'.format(date,grouping))
    #plt.show()
    plt.close('all')

def plot_polarplot(maxval, minval, highlatdata, midlatdata, val_name, keywrd, plot_title):
    fig = plt.figure()
    ax = fig.add_subplot(111, polar = True)

    cdict = {'red':  ((0.0, 0.8, 0.8),   # red at 0
                      (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                      (1.0, 0.0, 0.0)),  # no red at 1

            'green': ((0.0, 0.0, 0.0),   # no green at 0
                      (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                      (1.0, 0.0, 0.0)),  # green at 1

            'blue':  ((0.0, 0.0, 0.0),   # no blue at 0
                      (0.5, 1.0, 1.0),   # all channels set to 1.0 at 0.5 to create white
                      (1.0, 0.8, 0.8))   # no blue at 1
           }
    # Create the colormap
    GnRd = colors.LinearSegmentedColormap('GnRd', cdict)
    # Create the normalization (without norm only 0-1 values are allowed)
    norm = colors.Normalize(vmin=minval, vmax=maxval)

    # Plot each sector
    for i in range(4):
        highlatsign = ''
        midlatsign = ''

        # Plot mid lat
        color = GnRd(norm(midlatdata[i]))
        if i != 1:
            ax.bar(i*np.pi/2 - np.pi/2, 1, width=np.pi / 2, bottom=0,
                   color=color, edgecolor = color)
        else:
            ax.bar(i*np.pi/2 - np.pi/2, 1, width=np.pi / 2, bottom=0,
                   facecolor='#DEDEDE', hatch='/')
        if midlatdata[i] > 0:
            midlatsign = '+'
        else:
            midlatsign = '-'

        # Plot high lat
        color = GnRd(norm(highlatdata[i]))
        ax.bar(i*np.pi/2 - np.pi/2, 0.5, width=1 * np.pi / 2, bottom=0,
               color=color, edgecolor = color)
        if highlatdata[i] > 0:
            highlatsign = '+'
        else:
            highlatsign = '-'
            
        ax.annotate("{}{:.2f}".format(highlatsign, abs(highlatdata[i])), xy=[i*np.pi/2 - np.pi/2, 0.25], fontsize=14, ha='center',va='center')
        if i != 1:
            ax.annotate("{}{:.2f}".format(midlatsign, abs(midlatdata[i])), xy=[i*np.pi/2 - np.pi/2, 0.75], fontsize=14, ha='center',va='center')



        

    '''
    # Add in colorbar
    clb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=GnRd), ax=ax)
    clb.ax.set_title(r'{}'.format(val_name))
    '''
    ax.set_thetagrids((270, 0, 90, 180), labels = ['00 LT','6 LT', '12 LT', '18 LT'])
    ax.set_yticks([])
    #ax.set_yticklabels(['High Lat', '60$^\circ$', 'Mid Lat'])

    plt.plot([np.pi/4, np.pi/4], [0,1], 'k')
    plt.plot([3*np.pi/4, 3*np.pi/4], [0,1], 'k')
    plt.plot([5*np.pi/4, 5*np.pi/4], [0,1], 'k')
    plt.plot([7*np.pi/4, 7*np.pi/4], [0,1], 'k')
    plt.tick_params(labelsize=12)

    ax.spines['polar'].set_visible(False)
    angles = np.linspace(0,2*np.pi)
    r_vals = [0.5 for i in angles]
    plt.plot(angles,r_vals, 'k')
    r_vals = [1 for i in angles]
    plt.plot(angles,r_vals, 'k')

    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    fig.suptitle(plot_title)
    fig.savefig('{}_polarplot.png'.format(keywrd))
    plt.close()

def plot_per_smooth(thresholds, plot_type,dic):
    levels = ['60min', '30min', '1min', '15min']
    minval = 0
    maxval = 1
    if plot_type == 'heidke':
        maxval = 1
        minval = -1
        title = 'Heidke Skill Score for '
    elif plot_type == 'bias':
        minval = 0
        maxval=2
        title = 'Bias for '
    elif plot_type == 'delta':
        minval = -0.5
        maxval = 0.5
        title = 'Delta Skill for '

    for threshold in thresholds:
        for keywrd in levels:
            highlat_skill = []
            midlat_skill = []
            if plot_type == 'delta':
                if keywrd == '1min': continue
                else:
                    deltaskill_highlat = []
                    deltaskill_midlat = []
                    for i in range(5)[1:]:
                        delta_high = dic[threshold][keywrd]['highlat'][str(i)]['heidke'] - dic[threshold]['1min']['highlat'][str(i)]['heidke']
                        delta_mid =  dic[threshold][keywrd]['midlat'][str(i)]['heidke'] - dic[threshold]['1min']['midlat'][str(i)]['heidke']
                        highlat_skill += [delta_high]
                        midlat_skill += [delta_mid]
            else:
                for i in range(5)[1:]:
                    highlat_skill += [dic[threshold][keywrd]['highlat'][str(i)][plot_type]]
                    midlat_skill += [dic[threshold][keywrd]['midlat'][str(i)][plot_type]]
            
            plot_polarplot(maxval, minval, highlat_skill, midlat_skill, plot_type , '{0}_{1}'.format(keywrd,plot_type), title + keywrd)
          
def create_polarplot(thresholds, stations, starttimes):

    thresholds = [0.3]
    proc_dic, dic = create_dics(thresholds, stations, starttimes)
    
    midlat_stats = ['BEL', 'CLF', 'FMC', 'HAD', 'MEA', 'OTT', 'SIT', 'THY', 'WNG', 'DOU', 'FUR', 'HLP', 'PIN', 'STJ', 'UPS', 'BFE', 'ESK', 'GIM', 'NEW', 'PBQ', 'SUA', 'VAL', 'FCC', 'IRT', 'NGK', 'RAL', 'TAR', 'VIC']
    highlat_stats= ['ABK', 'BLC', 'BRW', 'BJN', 'CBB', 'CMO', 'DNB', 'DOB', 'EAG','FSP','SMI','HRN','IQA','STF','KEV','KUV','LER','LYR','NAQ','NAL','NRD','NUR','OUJ','THL','RAN','RES','SVS','TAL','AMK','TIK','YKC']

    plot_per_smooth(thresholds, 'heidke', dic)
    plot_per_smooth(thresholds, 'bias', dic)
    plot_per_smooth(thresholds, 'delta', dic)

    
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


def average_powerspectra(data, station, sector,smoothing,date):
    keywrds = ['1min','15min','30min','60min']
    fig = plt.figure()

    totals = []
    for i in range(4):
        curdata = np.array(data[i][1])
        keywrd = data[i][0]
        curdata = np.array(curdata)
        fourier_transform = np.fft.rfft(curdata)
        abs_fourier_transform = np.abs(fourier_transform)
        power_spectrum = np.square(abs_fourier_transform)

        time_step = 60
        frequency = np.fft.rfftfreq(curdata.size, time_step)
        idx = np.argsort(frequency)

        freq_sorted = frequency[idx]
        power_sorted = power_spectrum[idx]
        slow = freq_sorted < 0.001
        mid = np.logical_and(freq_sorted > 0.001,freq_sorted < 0.008)
        high = freq_sorted > 0.008
        slow_avg = np.average(power_sorted[slow])
        mid_avg = np.average(power_sorted[mid])
        high_avg = np.average(power_sorted[high])

        plt.loglog([0.001,0.004,0.008],[slow_avg,mid_avg,high_avg], label=keywrd)

        totals += [np.array([slow_avg,mid_avg,high_avg])]

    plt.legend()
    plt.ylabel(r'$[nT s]^2$')
    plt.xlabel('[Hz]')
    fig.savefig('plots/powerspectra/{}/combined/average_{}_{}_sector{}_{}_powerspectra.png'.format(date, date, station,sector,smoothing))

    totals[1] = totals[0] - totals[1]
    totals[2] = totals[0] - totals[2]
    totals[3] = totals[0] - totals[3]

    plt.close()

    return totals[1:]

def power_spectra(data, station, sector,smoothing,date):
    fig = plt.figure()

    return_power = []
    return_freq = []
    if len(data)==5:
        for i in range(5):
            curdata = np.array(data[i][1])
            keywrd = data[i][0]
            curdata = np.array(curdata)
            fourier_transform = np.fft.rfft(curdata)
            abs_fourier_transform = np.abs(fourier_transform)
            power_spectrum = np.square(abs_fourier_transform)
            time_step = 60
            frequency = np.fft.rfftfreq(curdata.size, time_step)
            idx = np.argsort(frequency)
            plt.loglog(frequency[idx], power_spectrum[idx],label = keywrd)
            return_power += [np.array(power_spectrum[idx])]
            return_freq += [np.array(frequency[idx])]
        plt.legend()

    else:
        data = np.array(data[1])
        fourier_transform = np.fft.rfft(data)
        abs_fourier_transform = np.abs(fourier_transform)
        power_spectrum = np.square(abs_fourier_transform)
        time_step = 60
        frequency = np.fft.rfftfreq(data.size, time_step)
        idx = np.argsort(frequency)
        plt.loglog(frequency[idx], power_spectrum[idx])

        return_power = [power_spectrum[idx]]
        return_freq = [frequency[idx]]
        #plt.show()

    plt.ylabel(r'$[nT s]^2$')
    plt.xlabel('[Hz]')

    if smoothing == 'combined':
        fig.savefig('plots/powerspectra/{}/combined/{}_{}_sector{}_{}_powerspectra.png'.format(date, date, station,sector,smoothing))
    else:
        fig.savefig('plots/powerspectra/{}/{}_{}_sector{}_{}_powerspectra.png'.format(date,date, station,sector,smoothing))
    plt.close()
    return return_freq, return_power



if __name__ == "__main__":
    #midlat_stats = ['BEL', 'CLF', 'FMC', 'HAD', 'MEA', 'OTT', 'SIT', 'THY', 'WNG', 'DOU', 'FUR', 'HLP', 'PIN', 'STJ', 'UPS', 'BFE', 'ESK', 'GIM', 'NEW', 'PBQ', 'SUA', 'VAL', 'FCC', 'IRT', 'NGK', 'RAL', 'TAR', 'VIC']
    midlat_stats= ['BEL','BOU','BFE','DOB','DOU','FRD','HAN','IRT','LER','NEW','NUR','OTT','SIT','STJ','UPS','VAL','VIC']
    #midlat_stats= ['BEL','BOU','BFE']
    #17

    #highlat_stats= ['ABK', 'BLC', 'BRW', 'BJN', 'CBB', 'CMO', 'DNB', 'DOB', 'EAG','FSP','SMI','HRN','IQA','STF','KEV','KUV','LER','LYR','NAQ','NAL','NRD','NUR','OUJ','THL','RAN','RES','SVS','TAL','AMK','TIK','YKC']
    #highlat_stats = ['ABK', 'YKC', 'IQA']
    highlat_stats = ['ABK','ATU','BJN','BET','DMH','DAW','EKP','IQA','HRN','LRV','MEA','NAQ','PBK','PBQ','PIN','THL','YKC','DOB']
    #highlat_stats = ['ABK','ATU','BJN']
    #18

    # this seciton does grouping by high lat and mid lat for single event
    #print('high lat')
    #grouping(outputdir, smdata, thresholds, highlat_stats, date, starttime)
    #print('midlat')
    #grouping(outputdir, smdata, thresholds, midlat_stats, date, starttime)



    # this section does grouping across all events for mid lat and high lat
    #starttimes = ['20061214120000','20010831000000','20050831100000','20100405000000']
    starttimes = ['20061214120000','20010831000000','20050831100000','20100405000000','20110805090000']
    #starttimes = ['20061214120000']


    stations = midlat_stats + highlat_stats
    thresholds = [0.3, 0.7, 1.1, 1.5] #nT/s
    '''
    print('high lat')
    cross_event_grouping(thresholds, highlat_stats, starttimes)
    print('mid lat')
    cross_event_grouping(thresholds, midlat_stats, starttimes)
    '''
    create_polarplot(thresholds, stations, starttimes)

    #main(sys.argv)