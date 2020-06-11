import numpy as np
import spacepy.time as spt


def magnetopause(swdata, theta, model='shue'):
    """Given input variables calculate magnetopause location

    Shue et al. model is axially symmetric in GSE
    """
    import spacepy.empiricals as emp
    #convert theta to local times
    localtime = (12 + theta/15) % 24
    # make sw input (P, Bz)
    try:
        fp = swdata['Flow_pressure'][0]
        bz = swdata['Bz_GSM'][0]
    except (TypeError, IndexError):
        fp = swdata['Flow_pressure']
        bz = swdata['Bz_GSM']
    indata = {'P': fp, 'Bz': bz}
    if model.lower() == 'shue':
        xyMP = emp.getMagnetopause(indata, LTs=localtime)[0]
        return xyMP[:, 0], xyMP[:, 1]


_chao_cf = {1: 11.1266, 2: 0.0010, 3: -0.0005,
              4: 2.5966, 5: 0.8182, 6: -0.017,
              7: -0.0122, 8: 1.3007, 9: -0.0049,
              10: -0.0328, 11: 6.047, 12: 1.029,
              13: 0.0231, 14:-0.002
              }


def _chao_alpha(swdata, cf):
    bzc = cf[13] if swdata['Bz_GSE'] >= 0 else cf[6]
    t1 = 1 + bzc*swdata['Bz_GSE']
    t2 = 1 + cf[7]*swdata['Flow_pressure']
    t3 = 1 + cf[10]*np.log(1 + swdata['Plasma_beta'])
    t4 = 1 + cf[14]*swdata['mach_ms']
    alp = cf[5] * t1 * t2 * t3 * t4
    return alp


def _chao_r0(swdata, cf):
    mms2 = swdata['mach_ms']**2
    bzc = cf[2] if swdata['Bz_GSE'] >= 0 else cf[3]
    part1 = cf[1] * (1 + bzc*swdata['Bz_GSE']) * (1 + cf[9]*swdata['Plasma_beta'])
    fracpart1 = (cf[8] - 1) * mms2 + 2
    fracpart2 = (cf[8] + 1) * mms2
    part2 = 1 + cf[4]*(fracpart1/fracpart2)*swdata['Flow_pressure']**(-1/cf[11])
    return part1 * part2


def bowshock(swdata, theta, model='chao2002'):
    """Given input variables calculate bow shock location

    Parameters
    ----------
    swdata : dict-like
        Input dictionary of solar wind parameters
        E.g., from spacepy.omni.get_omni(ticks, dbase='OMNI2hourly')
        Input dict must have:
        'Bx_GSE', 'By_GSE', 'Bz_GSE' [nT]
        'Ion_density' [cm^-3]
        'Plasma_temp' [K]
        'Plasma_bulk_speed' [km/s]
        'Plasma_beta'
        'Flow_pressure' [nPa]

    theta : float
        Angle (in degrees) from the nose of the bow shock

    Notes
    -----
    Chao et al. 2002 model is axially-symmetric in aberrated GSE coordinates
    """
    theta_r = np.deg2rad(theta)
    if model.lower() == 'chao2002':
        cf = _chao_cf
        swdata['mach_ms'] = mach_magnetosonic(swdata)
        # r_at_theta is radius of at given polar angle theta from
        # the aGSE (aberrated GSE) X-axis
        r_nought = _chao_r0(swdata, cf)
        flare = (1 + cf[12])/(1 + cf[12]*np.cos(theta_r))
        r_at_theta = r_nought*flare**_chao_alpha(swdata, cf)

    return r_at_theta


def get_BS_eq(swdata, angles=(-120, 120), npts=30):
    """Get x, y (GSE) for bow shock in the equatorial plane
    Input dict must have:
    'Bx_GSE', 'By_GSE', 'Bz_GSE' [nT]
    'Ion_density' [cm^-3]
    'Plasma_temp' [K]
    'Plasma_bulk_speed' [km/s]
    """
    angles = np.linspace(angles[0], angles[1], npts)
    angles_r = np.deg2rad(angles)
    radii = bowshock(swdata, angles)
    agsevec = np.zeros((npts, 3))
    agsevec[:, 0] = radii * np.cos(angles_r)
    agsevec[:, 1] = radii * np.sin(angles_r)
    gsevec = np.zeros_like(agsevec)
    for idx, row in enumerate(agsevec):
        gsevec[idx, ...] = agse_convert(row, swdata)
    return gsevec[:, 0], gsevec[:, 1]


def agse_convert(invec, swdata, insys='gse'):
    """Convert between GSE and aberrated GSE

    Using common definiton, e.g., Schwartz 1998
    (S.J. Schwartz, "Shock and Discontinuity Normals, Mach Numbers, and Related Parameters",
    In: Analysis Methods for Multi-Spacecraft Data, Eds.: G. Paschmann and P. Daly,
    ISSI Scientific Reports Series, ESA/ISSI, Vol. 1, ISBN 1608-280X, 1998, pp.249-270)

    Neglects transverse components of SW velocity
    """
    assert insys in ('gse', 'agse')
    alpha = np.arctan(30/swdata['Plasma_bulk_speed'])
    gse_to_agse = np.zeros((3,3), dtype=float)
    gse_to_agse[2, 2] = 1
    gse_to_agse[0, 0] = np.cos(alpha)
    gse_to_agse[1, 1] = np.cos(alpha)
    gse_to_agse[0, 1] = -np.sin(alpha)
    gse_to_agse[1, 0] = np.sin(alpha)

    if insys == 'gse':
        outvec = np.dot(gse_to_agse, invec)
    else:
        outvec = np.dot(gse_to_agse.T, invec)
    return outvec

def gse_to_agse(posvec, velvec):
    """conversion from GSE to aberrated GSE

    as given in Dmitriev et al., 2003

    posvec - 3-vector of position (GSE)
    velvec - 3-vector of velocity (GSE)
    """
    vx = velvec[0]
    vy = velvec[1]
    vz = velvec[2]
    x_gse = posvec[0]
    y_gse = posvec[1]
    z_gse = posvec[2]
    # rotate around Z_GSE by deltay
    deltay = np.arctan((vy + 30)/np.abs(vx))
    # rotate around aberrated Y by deltaz
    deltaz = np.arctan(vz/(np.sqrt(vx**2 + (vy+30)**2)))

    x_agse = (x_gse*np.cos(deltay) - y_gse*np.sin(deltay))*np.cos(deltaz) + z_gse*np.sin(deltay)
    y_agse = y_gse*np.cos(deltay) + x_gse*np.sin(deltay)
    z_agse = z_gse*np.cos(deltaz) + (x_gse*np.cos(deltay) - y_gse*np.sin(deltay))*np.sin(deltaz)

    outvec = np.array([x_agse, y_agse, z_agse])
    return outvec


def mach_magnetosonic(swdata):
    """
    Calculate the magnetosonic mach number

    Input dict must have:
    'Bx_GSE', 'By_GSE', 'Bz_GSE' [nT]
    'Ion_density' [cm^-3]
    'Plasma_temp' [K]
    'Plasma_bulk_speed' [km/s]
    """
    poly = 5/3  # Polytropic index

    # Angle between normal to bowshock and IMF vector
    theta_bn = np.pi/2
    # Field magnitude in nT
    Bmag = np.sqrt(swdata['Bx_GSE']**2 +
                   swdata['By_GSE']**2 +
                   swdata['Bz_GSE']**2)
    # Mass density. 1.16 is nominal accounting for 0.04 mass fraction of Helium
    nD = swdata['Ion_density']
    rho = 1.16*nD  # Mass density scaled to proton mass
    # Thermal pressure
    P_th = 1.6e-4**swdata['Plasma_temp']
    # Alfven speed
    v_a = 21.8*Bmag/np.sqrt(rho)
    # Sound speed
    v_s = np.sqrt(poly*P_th/rho)
    # Magnetosonic speed
    vavs2 = v_a**2 + v_s**2
    v_ms = (vavs2 + ((vavs2**2 - 4*v_a**2*v_s**2*np.cos(theta_bn)**2)**0.5)/2)**0.5
    # Then magnetosonic Mach number is the ratio of Vsw to Vms
    mach_ms = swdata['Plasma_bulk_speed']/v_ms
    return mach_ms


def plotMSBS(time=None, data=None, angles=(-125, 125), npts=30, mp=True, bs=True):
    """
    Plot bow shock and magnetopause

    Sample compressed time: 2012-10-02T03:00:00
    Sample "normal" time: 2012-07-22T12:00:00
    """
    import spacepy.omni as om
    import spacepy.empiricals as emp
    import spacepy.plot as splot
    import matplotlib.pyplot as plt
    if time is None and data is None:
        #use default time for "nominal" values
        time = spt.Ticktock('2012-10-02T03:00:00')
        omd = om.get_omni(time, dbase='OMNI2hourly')
    elif time is None and data is not None:
        omd = data
    elif time is not None and data is None:
        omd = om.get_omni(time, dbase='OMNI2hourly')
    else:
        raise ValueError('Use either "time" or "data" kwarg - not both.')
    fig, ax0 = plt.subplots(1)
    if mp:
        xmp, ymp = magnetopause(omd, np.linspace(angles[0], angles[1], npts))
        ax0.plot(xmp, ymp, label='Shue1997')
    if bs:
        x, y = get_BS_eq(omd, angles=angles, npts=npts)
        ax0.plot(x, y, label='Chao2002')
    ax0.set_aspect('equal')
    splot.dual_half_circle(ax=ax0)
    ax0.set_ylabel('Y$_{GSE}$ [R$_{E}$]')
    ax0.set_xlabel('X$_{GSE}$ [R$_{E}$]')
    ax0.legend()
    ax0.set_xlim([-40, 40])
    ax0.set_ylim([-45, 45])
    return fig, ax0

if __name__ == '__main__':
    import spacepy.empiricals as emp
    import matplotlib.pyplot as plt
    swdata89 = {}
    swdata89['Bx_GSE'] = 10
    swdata89['By_GSE'] = 15
    swdata89['Bz_GSE'] = -15
    swdata89['Bz_GSM'] = -15
    swdata89['Ion_density'] = 70
    swdata89['Plasma_bulk_speed'] = 700
    swdata89['Plasma_temp'] = emp.getExpectedSWTemp(swdata89['Plasma_bulk_speed'],
                                                    model='XB15', units='K')
    swdata89['Plasma_temp'] *= 1.2  # reduce by ~10% for CME-like temp 
    swdata89['Flow_pressure'] = swdata89['Plasma_bulk_speed']**2 *\
                                swdata89['Ion_density']*1.67621e-6
    bsq = swdata89['Bx_GSE']**2 + swdata89['By_GSE']**2 + swdata89['Bz_GSE']**2
    pm = 3.98e-4*bsq
    # Thermal pressure
    pth = 1.6e-4**swdata89['Plasma_temp']
    swdata89['Plasma_beta'] = pth/pm
    #print('beta = ' + '{}  ({}/{})'.format(swdata89['Plasma_beta'], pth, pm))

    fig, ax = plotMSBS(data=swdata89, mp=False, bs=True)
    # add IMP8 traj
    import os
    import datetime as dt
    import swtools
    import spacepy.datamanager as dman
    import spacepy.toolbox as tb
    datapath = os.path.abspath(os.path.join('..', 'ref_data', '1989'))
    # When plotting, use spacepy.datamanager to insert fill between contiguous regions
    mit_pl = swtools.readIMP8plasmafile(os.path.join(datapath, 'imp8.data.1989.060.090'))
    mit_t, mit_x = dman.insert_fill(np.asarray(mit_pl['time']), np.asarray(mit_pl['pos_gse'][:,0]))
    mit_t, mit_y = dman.insert_fill(np.asarray(mit_pl['time']), np.asarray(mit_pl['pos_gse'][:,1]))
    inds = tb.tOverlapHalf([dt.datetime(1989,3,12), dt.datetime(1989,3,15)], mit_t)
    ax.plot(mit_x[inds], mit_y[inds], label='MIT plasma')
    plt.show()
