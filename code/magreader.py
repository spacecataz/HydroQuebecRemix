'''
Derived class for reading CCMC magnetometer files into spacepy.pybats.bats.MagFile like class
    MagFile

See
https://ccmc.gsfc.nasa.gov/support/GEM_metrics_08/comments_ground.php
for details of the stations used, and
https://ccmc.gsfc.nasa.gov/RoR_WWW/pub/dBdt/out/
for the model output and observation data used for the
validation study
'''
import os
import sys
import dateutil.parser as dup
import numpy as np
import spacepy.pybats.bats


class CCMCMagFile(spacepy.pybats.bats.Mag):

    def __init__(self, fname, *args, **kwargs):
    #def __init__(self, nlines, time, gmvars=(), ievars=(), *args, **kwargs):
        from numpy import zeros

        super(spacepy.pybats.bats.Mag, self).__init__(*args, **kwargs)  # Init as PbData.

        with open(fname, 'r') as fh:
            contents = fh.readlines()

        header = [line.strip() for line in contents if not line[0].isdigit()]
        data = [line.strip() for line in contents if line[0].isdigit()]
        
        time = []
        nlines = len(data)
        self['dBn'] = np.zeros(nlines)
        self['dBe'] = np.zeros(nlines)
        self['dBd'] = np.zeros(nlines)
        for idx, line in enumerate(data):
            line = line.split()
            line_t = [int(el) for el in line[:5]]
            time.append(dup.parse('{0}-{1:02d}-{2:02d}T{3:02d}:{4:02d}'.format(int(line[0]),
                        int(line[1]), int(line[2]), int(line[3]), int(line[4]))))
            for key, col in [('dBn', 8), ('dBe', 9), ('dBd', 10)]:
                self[key][idx] = float(line[col])

        self['time'] = np.asarray(time)
        self.attrs['nlines'] = nlines
        iname = [line for line in header if line.startswith('#Inst')][0].split()[-1]
        self.attrs['instrument'] = iname

        self.calc_h()
        self.calc_dbdt()  # Ensure that same dB/dt method is used for data and sim


if __name__ == "__main__":
    # Run this when called as a script
    # We'll use it to demonstrate use.
    import matplotlib.pyplot as plt
    import spacepy.plot as splot

    indir = os.path.abspath(os.path.join('..', 'ref_data', 'ccmc'))
    fname = 'ykc_OBS_20050831.txt'  # Event 4, Yellowknife, from CCMC website

    # Instantiate reader
    data = CCMCMagFile(os.path.join(indir, fname))
    # plot delta-B and dB/dt
    col1 = 'royalblue'
    col2 = 'seagreen'
    splot.style('spacepy')
    fig, ax1 = plt.subplots(1, figsize=(10,4))
    ax2 = ax1.twinx()  # second axes object with shared x-axis
    ax2.grid(False)  # turn off the grid for the secondary axes
    data.add_plot('dBh', target=ax1, color=col1, ls='-')
    ax1.set_ylabel('dB$_h$ [nT]')
    ax1.set_ylim([0, 1000])
    data.add_plot('dBdth', target=ax2, color=col2, ls='-')
    ax2.set_ylabel('(dB/dt)$_h$ [nt/s]')
    ax2.set_ylim([0,3])
    ax2.yaxis.label.set_color(col2)  # Set color of second y-axis label
    ax2.tick_params(axis='y', colors=col2)  # Set color of 2nd y-axis ticks/ticklabels
    plt.tight_layout()  # fix layout so labels all fit...
    plt.show()
