import astropy.time as Time
from sunpy.time import parse_time
from stereo_spice.coordinates import StereoSpice
import scipy.io as sio
import pathlib
from scipy.io import readsav
import numpy as np
import tomo_huw as tomo

dir="/Users/hum2/data/stereo/secchi/pb/a/cor2/pythontest/"
files=sorted(pathlib.Path(dir).glob('*.dat'))

d={}
for file in files:
    v=readsav(file,verbose=False)
    tag_name=list(v.keys())
    values=list(v.values())
    print(tag_name[0])
    if tag_name[0]=="dates" or tag_name[0]=="files" or tag_name[0]=="filesorig":
        n=np.size(values[0])
        val = ["" for i in range(n)]
        for i in range(n):
            val[i]=values[0][i].decode("utf-8")
        values=[val]
    if tag_name[0]=="system" or tag_name[0]=="type" or tag_name[0]=="geometry":
        values=values[0].decode("utf-8")
    if np.size(values) > 1:
        values=values[0]
    d[tag_name[0]]=values

rmain=4.
sz=np.shape(d["im"])
nt=sz[0]
npa=sz[2]
nx=50
geom=tomo.tomo_make_geom(d,rmain,nt,npa,nx,dates={},spacecraft={})

# spice = StereoSpice()

# obstime=parse_time("2017/08/01 12:00:00")
# # obstime=Time("2017-08-01T12:00:00")
# sta=spice.get_lonlat(obstime,'sta','CARRINGTON',precess=False)
# print(sta)
# a=sta