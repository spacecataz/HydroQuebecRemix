# Tools for working with IMP8/ISEE data
import os
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
import spacepy.datamodel as dm
import spacepy.time as spt
import spacepy.datamanager as dman
import spacepy.plot as splot


def readIMP8plasmafile(fname):
    """
    ftp://space.mit.edu/pub/plasma/imp/fine_res/1989/


  yr  doy hh mm ss sc  decimal yr rg md    xse    yse     zse     ysm     zsm      speed     thermal speed     density      E/W angle     N/S angle
                                                                                mom  nonlin    mom  nonlin    mom nonlin   mom   best  thresh threshs
    """
    header = []
    with open(fname, 'r') as fh:
        while True:
            pos = fh.tell()
            line = fh.readline().strip()
            if not line:
                # empty line, skip
                continue
            if not line[0].isdigit():
                # this is header, save it
                header.append(line)
            else:
                # first line of data, roll back to start and pass to numpy
                fh.seek(pos)
                break
        data = np.loadtxt(fh)

    def toDate(yr, doy, hh, mm, ss):
        MM, DD = spt.doy2date(yr, doy)
        dates = dm.dmfilled(len(yr), fillval=None, dtype=object)
        for idx, (mon, day) in enumerate(zip(MM, DD)):
            dates[idx] = dt.datetime(yr[idx], mon, day, hh[idx], mm[idx], ss[idx])
        return dates

    region = data[:, 7]
    outdata = dm.SpaceData(attrs={'header': header, 'fname': fname})
    outdata['time'] = toDate(data[:, 0].astype(int), data[:, 1].astype(int),
                             data[:, 2].astype(int), data[:, 3].astype(int),
                             data[:, 4].astype(int))
    outdata['region'] = dm.dmarray(region)
    outdata['pos_gse'] = dm.dmarray(data[:, 9:12], attrs={'coord_sys': 'gse'})
    outdata['pos_gsm'] = dm.dmfilled(outdata['pos_gse'].shape, fillval=0,
                                     dtype=float, attrs={'coord_sys': 'gsm'})
    outdata['pos_gsm'][:, 0] = data[:, 9]
    outdata['pos_gsm'][:, 1:] = data[:, 12:14]
    outdata['speed'] = dm.dmarray(data[:, 14],
                                  attrs={'description': 'speed from moments'})
    # outdata['speed'][region > 2] = np.nan  # region 3 is sheath
    outdata['speed_nl'] = dm.dmarray(data[:, 15])
    vmask = outdata['speed_nl'] >= 9000
    # outdata['speed_nl'][region > 2] = np.nan  # region 3 is sheath
    outdata['speed_nl'][vmask] = np.nan  # region 3 is sheath
    outdata['n_dens'] = dm.dmarray(data[:, 18],
                                   attrs={'description': 'number density from moments'})
    outdata['n_dens_nl'] = dm.dmarray(data[:, 19])
    outdata['temp'] = 60.5*dm.dmarray(data[:, 16])**2
    outdata['temp_nl'] = 60.5*dm.dmarray(data[:, 17])**2
    outdata['data'] = data
    return outdata


def readIMPplasmaLANL(fname):
    """
    https://spdf.gsfc.nasa.gov/pub/data/imp/imp8/plasma_lanl/solarwind_2min/data/
    13 Flow Speed AP     E11.4      Relative Flow Speed    km/sec
       space             X1
    14 Flow Azimuth AP   E11.4      Relative Flow Angle    deg
       space             X1
    15 TRatio            E11.4      Alpha Temp/Proton Temp
       space             X1
    16 Temp              E11.4      Alpha Temp Anisotropy
       space             X1
    17 PSIT              E11.4      Alpha Pressure         deg
                                    Asymmetry Axis
    """
    with open(fname, 'r') as fh:
        # first line of data, roll back to start and read all
        data = fh.readlines()

    def toDate(ymd, hms):
        dates = dm.dmfilled(len(ymd), fillval=None, dtype=object)
        for idx, (p1, p2) in enumerate(zip(ymd, hms)):
            yr = p1//10000
            mon = (p1 - yr*10000)//100
            day = p1 - (yr*10000 + mon*100)
            hh = p2//10000
            mm = (p2 - hh*10000)//100
            ss = int(p2 - (hh*10000 + mm*100))
            dates[idx] = dt.datetime(int(yr)+1900, int(mon), int(day)) +\
                dt.timedelta(hours=int(hh), minutes=int(mm), seconds=ss)
        return dates
    data = np.array([line.strip().split()[1:] for line in data], dtype=float)
    outdata = dm.SpaceData()
    outdata['time'] = toDate(data[:, 0], data[:, 2])
    outdata['n_dens'] = dm.dmarray(data[:, 6])
    outdata['speed'] = dm.dmarray(data[:, 7])
    outdata['speed_rel'] = dm.dmarray(data[:, 13])
    outdata['alpha_proton_ratio'] = dm.dmarray(data[:, 12])
    return outdata


def readISEEmag(fname):
    """Read an ISEE magnetometer data file

    Data source: https://spdf.sci.gsfc.nasa.gov/pub/data/isee/isee3/magnetic_fields/1min_ascii_extracted/1min_hgi_1984_1990/
    """
    data = np.loadtxt(fname)

    def toDate(yr, doy, hh, mm):
        MM, DD = spt.doy2date(yr, doy)
        dates = dm.dmfilled(len(yr), fillval=None, dtype=object)
        for idx, (mon, day) in enumerate(zip(MM, DD)):
            dates[idx] = dt.datetime(yr[idx], mon, day, hh[idx], mm[idx])
        return dates
    outdata = dm.SpaceData()
    outdata['time'] = dm.dmarray(toDate(data[:, 0].astype(int), data[:, 1].astype(int),
                                        data[:, 2].astype(int), data[:, 3].astype(int)))
    outdata['B'] = dm.dmarray(data[:, 4:7], attrs={'coord_sys': 'GSE'})
    # replace bad values with NaN fill
    outdata['B'][outdata['B'] == 999.9] = np.nan
    return outdata


def plotEvent(st=dt.datetime(1989, 3, 11), en=dt.datetime(1989, 3, 15)):
    datapath = os.path.abspath(os.path.join('..', 'ref_data', '1989'))
    # When plotting, use spacepy.datamanager to insert fill between contiguous regions
    fig, ax = plt.subplots(5, sharex=True, figsize=(12, 8))
    # Do MIT plasma
    mit_pl = readIMP8plasmafile(os.path.join(datapath, 'imp8.data.1989.060.090'))
    mit_t, mit_sp = dman.insert_fill(np.asarray(mit_pl['time']), np.asarray(mit_pl['speed']))
    ax[0].plot(mit_t, mit_sp, label='MIT plasma')
    mit_t, mit_nd = dman.insert_fill(np.asarray(mit_pl['time']), np.asarray(mit_pl['n_dens']))
    ax[1].plot(mit_t, mit_nd, label='MIT plasma')

    # Do LANL plasma
    lanl_pl = readIMPplasmaLANL(os.path.join(datapath, '198903_imp8_lanl_sw_2min.asc'))
    lanl_t, lanl_sp = dman.insert_fill(np.asarray(lanl_pl['time']), np.asarray(lanl_pl['speed']))
    ax[0].plot(lanl_t, lanl_sp, label='LANL plasma')
    lanl_t, lanl_nd = dman.insert_fill(np.asarray(lanl_pl['time']), np.asarray(lanl_pl['n_dens']))
    ax[1].plot(lanl_t, lanl_nd, label='LANL plasma')

    # Do IMF from ISEE-3
    isee_ma = readISEEmag(os.path.join(datapath, '198903_isee3_mag03_1min.asc'))
    ax[2].plot(isee_ma['time'], isee_ma['B'][:, 0], label='ISEE-3 Bx')
    ax[2].plot(isee_ma['time'], isee_ma['B'][:, 1], label='ISEE-3 By')
    ax[2].plot(isee_ma['time'], isee_ma['B'][:, 2], label='ISEE-3 Bz')
    ax[2].axhline(linestyle='--', color=(0.3, 0.3, 0.3))

    # Add stuff from burton_test
    from burton import invert_example
    invert_example(axes=[ax[3], ax[4]], show=False)

    # Finalize
    ax[0].legend(loc='upper left', fancybox=True, framealpha=0.5)
    ax[0].set_ylabel('Speed\n[km/s]')
    ax[1].legend(loc='upper left', fancybox=True, framealpha=0.5)
    ax[1].set_ylabel('Number Density\n[cm$^{-3}$]')
    ax[2].legend(loc='upper left', fancybox=True, framealpha=0.5)
    ax[2].set_ylabel('ISEE3\nIMF [nT]')
    ax[2].set_ylim([-15, 15])
    splot.applySmartTimeTicks(ax[2], [st, en], dolabel=False)
    ax[2].set_xlim([st, en])
    plt.show()


def plotOrbits(st=dt.datetime(1989, 3, 11), en=dt.datetime(1989, 3, 15)):
    # Need to show orbits of IMP-8 and ISEE-3, along with estimated MP and bow shock locations...
    pass
