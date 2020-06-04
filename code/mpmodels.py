import numpy as np


def magnetopause():
    """Given input variables calculate magnetopause location
    """
    pass


def bowshock(swdata, model='chao2002'):
    """Given input variables calculate bow shock location
    """
    
    if model.lower() == 'chao2002':
        cf = {1: 11.1266, 2: 0.0010, 3: -0.0005,
              4: 2.5966, 5: 0.8182, 6: -0.017,
              7: -0.0122, 8: 1.3007, 9: -0.0049,
              10: -0.0328, 11: 6.047, 12: 1.029,
              13: 0.0231, 14:-0.002
              }
        def alpha(swdata, cf):
            bzc = cf[13] if swdata['Bz'] >= 0 else cf[6]
            t1 = 1 + bzc*swdata['Bz']
            t2 = 1 + cf[7]*swdata['Pdyn']
            t3 = 1 + cf[10]*np.log(1 + swdata['beta'])
            t4 = 1 + cf[14]*swdata['mach_ms']
            alp = cf[5] * t1 * t2 * t3 * t4
            return alp

        def r0(swdata, cf):
            bzc = cf[2] if swdata['Bz'] >= 0 else cf[3]
            part1 = cf[1] * (1 + cf[2]*swdata['Bz']) * (1 + cf[9]*swdata['beta'])
            fracpart = ((cf[8] - 1) * swdata['mach_ms']**2 + 2)/((cf[8] + 1) * swdata['mach_ms']**2)
            part2 = 1 + cf[4]*fracpart*swdata['Pdyn']**(-1/cf[11])
            return part1 * part2
        
        # r_at_theta is radius of at given polar angle theta from
        # the aGSE (aberrated GSE) X-axis
        r_at_theta = r0(swdata, cf) * ((1 + cf[12])/(1 + cf[12]*np.cos(theta)))
        r_at_theta = r_at_theta**alpha(swdata, cf)

    return r_at_theta


def gse_to_agse(posvec, velvec):
    """conversion from GSE to aberrated GSE

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
    poly = 5/3  # polytropic index

    # rho = mass density
    # B = magnetic field strength
    # P = solar wind pressure (n*k_B*(T_e+T_i))
    v_a = B/np.sqrt(rho*mu0)  # Alfven speed
    v_s = poly*P/rho  # Sound speed
    v_ms = np.sqrt(v_a**2 + v_s**2)
    mach_ms = v_sw/v_ms
    return mach_ms
