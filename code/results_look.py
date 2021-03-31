import spacepy
import spacepy.pybats.bats as bts
import sys
sys.path.append('/Users/sgraf/Desktop/SWMFtools')
sys.path.append('/Users/sgraf/Desktop/SWMFtools/dBdt')
import util
import matplotlib.pyplot as plt
import matplotlib
import supermag_parser

import spacepy.plot as splot
splot.style('spacepy')
plt.rcParams["legend.frameon"] = True
plt.rcParams["legend.facecolor"] = 'white'

def results_summary_update(log, geolog, show=True):
    """3-panel summary plot from log and geoindex files
    """
    fig, axes = plt.subplots(figsize=(10,10),nrows=3, ncols=1, sharex=True,gridspec_kw={'height_ratios': [2, 2, 1]})
    geolog.add_ae_quicklook(val='AU', plot_obs=True, target=axes[0], label='Unsmoothed',c='r')
    geolog.add_ae_quicklook(val='AL', plot_obs=True, target=axes[0],label=None, c='r')
    geolog.add_kp_quicklook(plot_obs=True, target=axes[1], label='Unsmoothed',c='r')
    log.add_dst_quicklook(plot_obs=True, target=axes[2], label='Unsmoothed',c='r')
    axes[0].set_xlabel('')
    axes[1].set_xlabel('')
    
    if show:
        plt.show()
    return fig, axes

dates = ['20061214','20010830','20050831','20100405','20110805','20150316']
for date in dates:
    hour_logs = util.load_logs('/Users/sgraf/Desktop/SWMF_analysis/outputs/{}/hour/'.format(date))
    orig_logs = util.load_logs('/Users/sgraf/Desktop/SWMF_analysis/outputs/{}/unsmoothed/'.format(date))
    thirty_logs = util.load_logs('/Users/sgraf/Desktop/SWMF_analysis/outputs/{}/30min/'.format(date))
    thirty_geo = util.load_logs('/Users/sgraf/Desktop/SWMF_analysis/outputs/{}/30min'.format(date),logtype='geo')
    hour_geo = util.load_logs('/Users/sgraf/Desktop/SWMF_analysis/outputs/{}/hour/'.format(date),logtype='geo')
    orig_geo = util.load_logs('/Users/sgraf/Desktop/SWMF_analysis/outputs/{}/unsmoothed/'.format(date),logtype='geo')

    fig, axes = results_summary_update(orig_logs,orig_geo,show=False)

    hour_geo.add_ae_quicklook(val='AU', plot_obs=False, target=axes[0], label='Hourly', c='b')
    hour_geo.add_ae_quicklook(val='AL', plot_obs=False, target=axes[0],label='Test', c='b')
    hour_geo.add_kp_quicklook(plot_obs=False, target=axes[1], label='Hourly',c='b')
    hour_logs.add_dst_quicklook(plot_obs=False, target=axes[2], label='Hourly',c='b')

    thirty_geo.add_ae_quicklook(val='AU', plot_obs=False, target=axes[0], label='30min',c='g')
    thirty_geo.add_ae_quicklook(val='AL', plot_obs=False, target=axes[0],label=None,c='g')
    thirty_geo.add_kp_quicklook(plot_obs=False, target=axes[1], label='30min',c='g')
    thirty_logs.add_dst_quicklook(plot_obs=False, target=axes[2], label='30min',c='g')

    axes[0].get_legend().remove()
    axes[2].get_legend().remove()

    axes[0].set_ylabel('AU/AL (nT)')

    fig.suptitle('Summary Plot')
    fig.subplots_adjust(top=0.88)
    #plt.show()
    plt.savefig('{}_au_al_summary_plot_overlay.png'.format(date))
    #fig, axes = results_summary_update(hour_logs,hour_geo,show=False)
    #fig.suptitle('Hourly Smoothed Summary Plot')
    #fig.subplots_adjust(top=0.88)
    #plt.savefig('{}_au_al_summary_plot_hourly.png'.format(date))
    #fig, axes = results_summary_update(thirty_logs,thirty_geo,show=False)
    #fig.suptitle('30min Smoothed Summary Plot')
    #fig.subplots_adjust(top=0.88)
    #plt.savefig('{}_au_al_summary_plot_30min.png'.format(date))
'''
diff_dst = orig_logs['dst_sm'] - hour_logs['dst_sm']
fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)
orig_logs.add_dst_quicklook(plot_obs=True, target=axes[0])
hour_logs.add_dst_quicklook(plot_obs=True, target=axes[1])
axes[2].plot(orig_logs['time'],diff_dst)
axes[0].set_title('Dst of Unsmoothed Data')
axes[1].set_title('Dst of Hourly Smoothed Data')
axes[2].set_title('Difference between unsmoothed and hourly Dst')
fig.tight_layout()
'''
'''

fig,ax = plt.subplots()
ax.plot(orig_logs.obs_dst['time'], orig_logs.obs_dst['dst'], '--',color='black',label='Observed DST')
ax.plot(orig_logs['time'], orig_logs['dst_sm'],label='Unsmoothed DST')
ax.plot(hour_logs['time'], hour_logs['dst_sm'],label='Hourly DST')
plt.xlim(hour_logs['time'][0], hour_logs['time'][-1])
leg = ax.legend()
'''