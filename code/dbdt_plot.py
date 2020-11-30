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
import glob
import spacepy.plot as splot
splot.style('spacepy')

timerange = [dt.datetime(2015,3,16,16,00,00), dt.datetime(2015,3,19,00,00,00)]

def main(args):
    starttime = args[1]
    date = starttime[:8]
    year = date[:4]
    month = date[4:6]
    day = date[6:8]
    hour = starttime[8:10]
    minute = starttime[10:12]
    second = starttime[12:]

    endtime = args[2]
    starttime = dt.datetime.strptime(starttime, '%Y%m%d%H%M%S')

    global timerange
    timerange = [starttime, dt.datetime.strptime(endtime, '%Y%m%d%H%M%S')]
    print(timerange)
    outputdir = './outputs/{}/'.format(date)

    smdata = supermag_parser.supermag_parser('./supermag_data/{0}{1}{2}-supermag.txt'.format(year,month,day))
    stations = ['YKC', 'MEA', 'NEW', 'FRN', 'IQA', 'PBQ', 'OTT', 'FRD', 'VIC']
    make_plots(smdata,stations,outputdir,date, starttime)

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

def make_plots(smdata, stations, outputdir, date, starttime):
    unsmoothedfiles = glob.glob(outputdir + '/unsmoothed/mag*mag')
    hourlyfiles = glob.glob(outputdir + 'hour/mag*mag')
    files_30 = glob.glob(outputdir + '30min/mag*mag')

    #unsmoothedfiles = ['magnetometers_e20150316-160000.mag','magnetometers_e20150317-160000.mag','magnetometers_e20150318-160000.mag']
    #hourlyfiles = ['magnetometers_e20150316-160000_hour.mag','magnetometers_e20150317-160000_hour.mag']

    for stat in stations:
        fig, (ax1, ax2, ax3) = plt.subplots(figsize=(10,9),nrows=3, sharey=True)
        fig2, (b1, b2, b3) = plt.subplots(figsize=(10,9),nrows=3, sharey=True)

        count = 0
        for fname in unsmoothedfiles:
            origdata = spacepy.pybats.bats.MagFile(fname)
            origdata.calc_h()
            origdata.calc_dbdt()
            #if stat not in smdata.station: continue
            subset = origdata[stat]
            simtime = subset['time'][::6]
            dBdte = decimate(subset['dBdte'], 6)
            dBdtn = decimate(subset['dBdtn'], 6)
            dBdth = np.array([linalg.norm([dBdtn[i], dBdte[i]]) for i in range(len(dBdtn))])

            dBe = decimate(subset['dBe'], 6)
            dBn = decimate(subset['dBn'], 6)
            Bh = np.array([linalg.norm([dBn[i], dBe[i]]) for i in range(len(dBn))])

            b1.plot(simtime, Bh, 'b-', alpha=0.4)
            ax1.plot(simtime, dBdth, 'b-', alpha=0.4)

            run20, t20 = tb.windowMean(dBdth, time=simtime, winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
            brun20, bt20 = tb.windowMean(Bh, time=simtime, winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
            if count == 0:
                ax1.plot(t20, run20, marker='o', color='b', linestyle='none', markersize=3, label='SWMF Unsmoothed')
                b1.plot(bt20, brun20, marker='o', color='b', linestyle='none', markersize=3, label='SWMF Unsmoothed')
                count = 1
            else:
                ax1.plot(t20, run20, marker='o', color='b', linestyle='none', markersize=3)
                b1.plot(bt20, brun20, marker='o', color='b', linestyle='none', markersize=3)
        count = 0

        if len(hourlyfiles) > 0:
            for fname in hourlyfiles:
                origdata = spacepy.pybats.bats.MagFile(fname)
                origdata.calc_h()
                origdata.calc_dbdt()

                subset = origdata[stat]
                simtime = subset['time'][::6]
                dBdte = decimate(subset['dBdte'], 6)
                dBdtn = decimate(subset['dBdtn'], 6)
                dBdth = np.array([linalg.norm([dBdtn[i], dBdte[i]]) for i in range(len(dBdtn))])

                dBe = decimate(subset['dBe'], 6)
                dBn = decimate(subset['dBn'], 6)
                Bh = np.array([linalg.norm([dBn[i], dBe[i]]) for i in range(len(dBn))])


                ax3.plot(simtime, dBdth, 'b-', alpha=0.4)
                b3.plot(simtime, Bh, 'b-', alpha=0.4)

                run20, t20 = tb.windowMean(dBdth, time=simtime, winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
                brun20, bt20 = tb.windowMean(Bh, time=simtime, winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
                if count == 0:
                    ax3.plot(t20, run20, marker='o', color='b', linestyle='none', markersize=3, label='SWMF Hourly')
                    b3.plot(bt20, brun20, marker='o', color='b', linestyle='none', markersize=3, label='SWMF Hourly')
                    
                    count = 1
                else:
                    ax3.plot(t20, run20, marker='o', color='b', linestyle='none', markersize=3)
                    b3.plot(bt20, brun20, marker='o', color='b', linestyle='none', markersize=3)
        count = 0
        if len(files_30) > 0:
            for fname in files_30:
                origdata = spacepy.pybats.bats.MagFile(fname)
                origdata.calc_h()
                origdata.calc_dbdt()

                subset = origdata[stat]
                simtime = subset['time'][::6]
                dBdte = decimate(subset['dBdte'], 6)
                dBdtn = decimate(subset['dBdtn'], 6)
                dBdth = np.array([linalg.norm([dBdtn[i], dBdte[i]]) for i in range(len(dBdtn))])

                dBe = decimate(subset['dBe'], 6)
                dBn = decimate(subset['dBn'], 6)
                Bh = np.array([linalg.norm([dBn[i], dBe[i]]) for i in range(len(dBn))])


                ax2.plot(simtime, dBdth, 'b-', alpha=0.4)
                b2.plot(simtime, Bh, 'b-', alpha=0.4)

                run20, t20 = tb.windowMean(dBdth, time=simtime, winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
                brun20, bt20 = tb.windowMean(Bh, time=simtime, winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
                if count == 0:
                    ax2.plot(t20, run20, marker='o', color='b', linestyle='none', markersize=3, label='SWMF 30min')
                    b2.plot(bt20, brun20, marker='o', color='b', linestyle='none', markersize=3, label='SWMF 30min')
                    
                    count = 1
                else:
                    ax2.plot(t20, run20, marker='o', color='b', linestyle='none', markersize=3)
                    b2.plot(bt20, brun20, marker='o', color='b', linestyle='none', markersize=3)

        # add in observed
        if stat not in smdata.station: continue

        else:
            smstat = smdata.station[stat]

            Bdoth = np.array([linalg.norm(smstat['Bdot'][i,:2]) for i in range(len(smstat['Bdot']))])

            Bh = np.array([linalg.norm(smstat['B'][i,:2]) for i in range(len(smstat['B']))])

            ax1.plot(smstat['time'], Bdoth, 'r-', alpha=0.4)  
            ax2.plot(smstat['time'], Bdoth, 'r-', alpha=0.4)
            ax3.plot(smstat['time'], Bdoth, 'r-', alpha=0.4)

            b1.plot(smstat['time'], Bh, 'r-', alpha=0.4)  
            b2.plot(smstat['time'], Bh, 'r-', alpha=0.4)
            b3.plot(smstat['time'], Bh, 'r-', alpha=0.4)

            obs20, t20 = tb.windowMean(Bdoth, time=smstat['time'], winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
            bobs20, bt20 = tb.windowMean(Bh, time=smstat['time'], winsize=dt.timedelta(minutes=20), overlap=dt.timedelta(0), st_time=starttime, op=np.max)
            
            # Plot on dBdt plots
            ax1.plot(t20, obs20, marker='x', color='r', linestyle='none', markersize=3, label='Obs')
            ax2.plot(t20, obs20, marker='x', color='r', linestyle='none', markersize=3, label='Obs')
            ax3.plot(t20, obs20, marker='x', color='r', linestyle='none', markersize=3, label='Obs')

            # Plot on B plots
            b1.plot(bt20, bobs20, marker='x', color='r', linestyle='none', markersize=3, label='Obs')
            b2.plot(bt20, bobs20, marker='x', color='r', linestyle='none', markersize=3, label='Obs')
            b3.plot(bt20, bobs20, marker='x', color='r', linestyle='none', markersize=3, label='Obs')

            
        adjust_plots(ax1, '1-min dB/dt$_{H}$ [nT/s]', Zero=False, xlab=True)
        ax1.set_xlabel('Time')
        ax1.set_title('Unsmoothed')
        adjust_plots(ax2, '1-min dB/dt$_{H}$ [nT/s]', Zero=False, xlab=True)
        ax2.set_xlabel('Time')
        ax3.set_title('Hourly Smoothed')
        adjust_plots(ax3, '1-min dB/dt$_{H}$ [nT/s]', Zero=False, xlab=True)
        ax3.set_xlabel('Time')
        ax2.set_title('30min Smoothed')
        ax1.legend()
        ax2.legend()
        ax3.legend()

        adjust_plots(b1, '1-min B$_{H}$ [nT]', Zero=False, xlab=True)
        b1.set_xlabel('Time')
        b1.set_title('Unsmoothed')
        adjust_plots(b2, '1-min B$_{H}$ [nT]', Zero=False, xlab=True)
        b2.set_xlabel('Time')
        b3.set_title('Hourly Smoothed')
        adjust_plots(b3, '1-min B$_{H}$ [nT]', Zero=False, xlab=True)
        b3.set_xlabel('Time')
        b2.set_title('30min Smoothed')
        b1.legend()
        b2.legend()
        b3.legend()

        #splot.applySmartTimeTicks(ax, subset['time'], dolimit=True)

        fig.suptitle(stat)
        fig.tight_layout()
        fig.subplots_adjust(top=0.88)

        fig2.suptitle(stat)
        fig2.tight_layout()
        fig2.subplots_adjust(top=0.88)
        #plt.show()
        fig.savefig('/Users/sgraf/Desktop/SWMF_analysis/plots/dBdt/{0}_dBdt_comp_obs_{1}.png'.format(date, stat))
        fig2.savefig('/Users/sgraf/Desktop/SWMF_analysis/plots/dBdt/{0}_B_comp_obs_{1}.png'.format(date, stat))
        #plt.close()
if __name__ == "__main__":
    main(sys.argv)