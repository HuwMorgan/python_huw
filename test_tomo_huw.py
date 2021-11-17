import astropy.time as Time
from sunpy.time import parse_time
from stereo_spice.coordinates import StereoSpice
import scipy.io as sio
import pathlib
from scipy.io import readsav
import numpy as np
import tomo_huw as tomo
import json
import pickle
import matplotlib.pyplot as plt
import general_huw as gen
import coordinates_huw as coord
import constants_huw as const
import math
import functions_huw as func

# lon=coord.make_coordinates(360,[-1,1])*math.pi
# lat=0.1
# l=7
# m=6
# sphnow=func.make_spher_harm(7,lon,lat,userlm=[l,m])


dir="/Users/hum2/data/stereo/secchi/pb/a/cor2/pythontest/"
files=sorted(pathlib.Path(dir).glob('*.dat'))
savename=dir+'geom.pkl'
savenamesphrecon=dir+'sphrecon.pkl'
savenamesphdata=dir+'sphdata.pkl'

afile=open(savename,"rb")
geom=pickle.load(afile)
afile.close()

sphrecon,sphdata=tomo.tomo_prep_sph(geom)

afile=open(savenamesphrecon,"wb")
pickle.dump(sphrecon,afile)
afile.close()

afile=open(savenamesphdata,"wb")
pickle.dump(sphdata,afile)
afile.close()

a=b

# FOLLOWING TO OPEN TOMO DATACUBE STRUCTURE (IDL SAVE FILE) AND RUN TOMO_MAKE_GEOM, SAVING TO GEOM.PKL IN DIRECTORY ABOVE
# d={}
# for file in files:
#     v=readsav(file,verbose=False)
#     tag_name=list(v.keys())
#     values=list(v.values())
#     print(tag_name[0])
#     if tag_name[0]=="dates" or tag_name[0]=="files" or tag_name[0]=="filesorig":
#         n=np.size(values[0])
#         val = ["" for i in range(n)]
#         for i in range(n):
#             val[i]=values[0][i].decode("utf-8")
#         values=[val]
#     if tag_name[0]=="system" or tag_name[0]=="type" or tag_name[0]=="geometry":
#         values=values[0].decode("utf-8")
#     if np.size(values) > 1:
#         values=values[0]
#     d[tag_name[0]]=values

# rmain=4.
# sz=np.shape(d["im"])
# nt=sz[0]
# npa=sz[2]
# nx=50
# geom=tomo.tomo_make_geom(d,rmain,nt,npa,nx,dates={},spacecraft={})

# afile=open(savename,"wb")
# pickle.dump(geom,afile)
# afile.close()
