import os
import pickle
import numpy as np
from rebin import rebin
from pathlib import Path

# from tomo.time_huw import anytim2cal
# from tomo.general_huw import sizearr
# from tomo.tomo_prep_sph import tomo_prep_sph
# from tomo.tikhonov_opt_tomo import tikhonov_opt_tomo
# from tomo.timerange2datacube import timerange2datacube
# from tomo.tikhonov_search_tomo import tikhonov_search_tomo
# from tomo.tomo_prep_data import tomo_prep_data
# from tomo.tomo_huw import tomo_make_geom

from time_huw import anytim2cal
from general_huw import sizearr
from tomo_huw import tomo_prep_sph
from tomo_huw import tikhonov_opt_tomo
from tomo_huw import timerange2datacube
from tomo_huw import tikhonov_search_tomo
from tomo_huw import tomo_prep_data
from tomo_huw import tomo_make_geom

def tomo_data_2021(rmain, middate, overwrite=False):
    """ Demonstrates tomography processing and
        writes data to files in Python's Pickle format.
    
        Usage
        -----
        >>> tomo_data_2021(4., '2017/07/07')
        Parameters
        ----------
        `rmain`: `float`
            distance from Sun
        `middate`: `str`
            date for tomography
        `overwrite`: `bool`, optional, default `False`
            whether to perform tomography process if the file already exists
        Returns
        -------
            None
        Note
        ----
        Simplified version of tomo_data_2019, Huw 2021/11/24.
    """

    SECS_IN_DAY = 86400
    BASE_DIR = os.path.expanduser('~')

    dirdate = anytim2cal(middate, form=8, date_only=True)  # format date in yyyymmdd form e.g. 20170707, see
    dir_ = os.path.join(BASE_DIR, 'data', 'CorTom',
                        'tomography', 'cor2a', dirdate[0])  # main save directory
    if not Path(dir_).exists():
        os.makedirs(dir_,exist_ok=True)

    stereo = 'a'
    instr = 'cor2a'

    htstr = (str(round(10*rmain)/10.).replace(' ', ''))[:4]  # strmid(0,4) -- string of rmain, length 4 e.g. "4.50"

    savename = os.path.join(dir_, '_'.join(['tomo', htstr, dirdate[0], 'tomo.pkl']))  # name of file to be saved
    if not overwrite and os.path.isfile(savename):
        print('File exists', os.path.basename(savename))
        return

    # some fixed parameters
    mintrange = 10 * SECS_IN_DAY
    maxtimegap = 3 * SECS_IN_DAY
    norder = 8#22
    nlon = 180#540
    nlat = 90#270
    shrinkfactusr = 2
    nl = 10#25
    ndens = 5#20
    large=True

    print('Preparing datacube for date', middate)
    d = timerange2datacube(middate=middate,
                           mintrange=mintrange, maxtimegap=maxtimegap)
    if type(d) is not dict:  # timerange2datacube returns null
        raise TypeError('No data available (tomo_data_2021)')

    print('Extracting data at height', rmain, 'for date', middate)
    d= tomo_prep_data(d, rmain=rmain,shrinkfactusr=shrinkfactusr)
    files = d['files']
    filesorig = d['filesorig']

    print('Calculating geometry values at height', rmain, 'for date', middate)
    geom = tomo_make_geom(d,large=large)

    print('Preparing spherical harmonical basis for date', middate)
    sphrecon, sphdata = tomo_prep_sph(geom, nlon, nlat, norder=norder)
    lon = sphrecon['lon']
    lat = sphrecon['lat']

    print('Finding optimal reconstruction parameters for date', middate)
    
    t = tikhonov_search_tomo(d, sphrecon, sphdata,geom, nl=nl, ndens=ndens)  

    print('Applying tomography with optimal parameters for date', middate)
    c = tikhonov_opt_tomo(t)

    print('Cleaning observations and repeating tomography...')
    npa, nt, nsph = sizearr(sphdata['sph']) 
    bmod = np.sum(np.broadcast_to(np.expand_dims(c, axis=(0,1)), (npa, nt, nsph)) * sphdata['sph'], 2)
    bobs = d["im"]
    relnoise = d["relnoise"]
    dev = (bobs-bmod) / bobs
    thresh = 2 * np.nanstd(dev)
    index = abs(dev) >= thresh
    count = np.count_nonzero(index)
    if count > 0:
        print(count, 'bad pixels found')
        bobs[index] = np.nan
        relnoise[index] = np.nan
        d["im"]=bobs
        d["relnoise"]=relnoise
        t = tikhonov_search_tomo(d, sphrecon, sphdata,geom, nl=nl, ndens=ndens)
        c = tikhonov_opt_tomo(t)
        bmod = np.sum(np.broadcast_to(np.expand_dims(c, axis=(0,1)), (npa, nt, nsph)) * sphdata['sph'], 2)
    else:
        print('No outlying pixels found')

    nlon, nlat, nsph = sizearr(sphrecon['sph']) 
    dens = np.sum(np.broadcast_to(np.expand_dims(c, axis=(0,1)), (nlon, nlat, nsph)) * sphrecon['sph'], 2)
    

    d = {'rmain': rmain, 'date': middate, 'bobs': bobs,'bmod':bmod,
         'relnoise': relnoise, 'shrinkfactusr': shrinkfactusr,
         'instr': instr, 'norder': norder, 'nlon': nlon, 'nlat': nlat,
         'lon':lon, 'lat':lat, 't': t, 'c': c, 'dens': dens,
         'files': files, 'filesorig': filesorig}

    with open(savename, 'wb') as f:
        pickle.dump(d, f)
    
    return d