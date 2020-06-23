import datetime as dt
import numpy as np
import spacepy.datamodel as dm
import spacepy.toolbox as tb
from spacepy.pybats import kyoto as kyo
import spacepy.empiricals as emp

from burton import inverseOBrienMcPherron


def angle_between_vectors(vec1, vec2):
    """Angle between two vectors

    E.g., angle between x-axis and a position vector [3,4,5]:
    >>> vec1 = np.array([1, 0, 0])
    >>> vec2 = np.array([3, 4, 5])
    >>> angle_between_vectors(vec1, vec2)
    """
    vec1_norm = vec1/np.linalg.norm(vec1)
    vec2_norm = vec2/np.linalg.norm(vec2)
    dot_v1_v2 = np.dot(vec1_norm, vec2_norm)
    angle_deg = np.rad2deg(np.arccos(dot_v1_v2))
    return angle_deg


def getKondrashovSW(fn='../ref_data/1989/Kondrashov_SSA_19890312_19890314.txt'):
    data = dm.readJSONheadedASCII(fn, convert=True)
    return data


def rectify(invec):
    """Rectify input series so that negative values are zero"""
    invec[invec <= 0] = 0
    return invec


def calc_P(indata, v_var='V_sw', n_var='Den_P'):
    """Calculate Pdyn using V and n"""
    indata['Pdyn'] = indata[v_var]**2 * indata[n_var] * 1.67621e-6
    return indata


def downscaleSW(indata):
    """Downscale from hourly to 1-minute solar wind"""
    pass


def initialSW():
    """Construct initial estimate of hourly solar wind parameters"""
    st = dt.datetime(1989, 3, 12)
    en = dt.datetime(1989, 3, 15)
    hourly = getKondrashovSW()
    # Keep IMF By and n, set V to Boteler/Nagatsuma
    t_ssc1 = dt.datetime(1989, 3, 13, 1, 27)
    t_ssc2 = dt.datetime(1989, 3, 13, 7, 43)
    t_turn = dt.datetime(1989, 3, 14, 2)
    # Set Bx positive, in accordance with ISEE-3 data
    hourly['Bx'] = dm.dmfilled(hourly['By'].shape, fillval=3)
    # Before first SSC
    inds_before_1 = tb.tOverlapHalf([dt.datetime(1989, 3, 12), t_ssc1], hourly['DateTime'])
    hourly['V_sw'][inds_before_1] = 400
    # Between first and ssecond SSC
    inds_between_12 = tb.tOverlapHalf([t_ssc1, t_ssc2], hourly['DateTime'])
    hourly['V_sw'][inds_between_12] = 550
    # IMF turns north around 1989-03-14T02:00:00 according to inverse Burton and Kondrashov
    inds_mainphase = tb.tOverlapHalf([t_ssc2, t_turn], hourly['DateTime'])
    hourly['V_sw'][inds_mainphase] = 983
    # Then have speed decay towards IMP-8 measurement which is ballpark 820 km/s
    inds_rest = tb.tOverlapHalf([t_turn, hourly['DateTime'][-1]], hourly['DateTime'])
    hourly['V_sw'][inds_rest] = np.linspace(983, 820, len(inds_rest))
    # Get "Kondrashov VBs" using V from Boteler/Nagatsuma
    hourly['VBs_K'] = 1e-3 * hourly['V_sw'] * rectify(-1*hourly['Bz'])  # mV/m
    # Now get VBs from inverse Burton
    ky_dat = kyo.fetch('dst', (st.year, st.month, st.day), (en.year, en.month, en.day))
    inds = tb.tOverlapHalf([st, en], ky_dat['time'])
    # Pressure correct here using n and V
    hourly = calc_P(hourly)
    hourly['Dst'] = ky_dat['dst'][inds]
    dst_star = emp.getDststar(hourly, model='OBrien')
    hourly['VBs_OB'] = 1e-3*inverseOBrienMcPherron(dst_star)
    # Make new Bz from VBs_OB
    hourly['Bz_OB'] = -1e3 * hourly['VBs_OB']/hourly['V_sw']  # nT

    return hourly
