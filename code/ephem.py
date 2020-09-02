import os
import glob
import datetime as dt
from functools import partial
import numpy as np
import dateutil.parser as dup
from spacepy.pybats import SatOrbit
import spacepy.coordinates as spc

datadir = os.path.abspath('../ref_data/Ephem/1989')

for sat in ['GOES6', 'GOES7']:
    fns = glob.glob(os.path.join(datadir, '*{}*').format(sat))
    fns = sorted(fns)
    # read and concatenate position arrays
    with open(fns[0], 'r') as fh:
        data = fh.readlines()
    for fn in fns[1:]:
        with open(fn, 'r') as fh:
            data.extend(fh.readlines())
    # drop header lines
    data = [line.strip().split() for line in data if not line.startswith('#')]

    # grab times and positions
    parsetime = partial(dup.parse, ignoretz=True)
    tmp = np.asarray(data)
    time = tmp[:, 0]
    pos = tmp[:, 1:].astype(float)
    time = np.vectorize(parsetime)(time)
    rll_order = [2, 0, 1]  # reorder columns for coords
    geo_sph = spc.Coords(pos[:, rll_order], 'SPH', 'sph')
    geo_car = geo_sph.convert('SPH', 'car')

    # make and fill SatOrbit for writing
    satdata = SatOrbit()
    satdata['time'] = time
    satdata['xyz'] = geo_car.data.T  #SatOrbit is 3xN, not Nx3
    satdata.attrs['coor'] = 'GEO'
    outfname = os.path.join(os.path.abspath('..'), 'swmf_input', 'Event19890312',
                            '{}.dat'.format(sat))
    satdata.attrs['file'] = outfname
    satdata.write()
