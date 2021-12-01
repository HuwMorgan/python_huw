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
import os

# lon=coord.make_coordinates(360,[-1,1])*math.pi
# lat=0.1
# l=7
# m=6
# sphnow=func.make_spher_harm(7,lon,lat,userlm=[l,m])

rmain=4.
dir="/Users/hum2/data/sweep/"
overwrite_data=False
overwrite_geom=False
overwrite_sphrecon=False
overwrite_fit=True
overwrite_coeff=True

savenamed=dir+'d.pkl'
savenamegeom=dir+'geom.pkl'
savenamesphrecon=dir+'sphrecon.pkl'
savenamesphdata=dir+'sphdata.pkl'
savefit=dir+"tikhonov_search.pkl"
savecoeff=dir+"tikhonov_opt.pkl"

# FOLLOWING TO OPEN DATA DATACUBE STRUCTURE (IDL SAVE FILE) AND RUN TOMO_MAKE_GEOM, SAVING TO GEOM.PKL IN DIRECTORY ABOVE
if os.path.exists(savenamed) and not overwrite_data:
    print("Reading ",savenamed)
    afile=open(savenamed,"rb")
    d=pickle.load(afile)
    afile.close()
else:
    v=readsav(dir+'d.dat',verbose=False)
    tag_name=list(v.keys())
    values=list(v.values())
    d={}
    for i in np.arange(np.size(tag_name)):
        print(tag_name[i])
        if tag_name[i]=="dates" or tag_name[i]=="files" or tag_name[i]=="filesorig":
            n=np.size(values[i])
            val=[values[i][j].decode("utf-8") for j in range(n)]
        elif tag_name[i]=="system" or tag_name[i]=="type" or tag_name[i]=="geometry":
            val=values[i].decode("utf-8")
        else:
            val=values[i]
        d[tag_name[i]]=val
    afile=open(savenamed,"wb")
    pickle.dump(d,afile)
    afile.close()

# FOLLOWING TO OPEN DATA DATACUBE STRUCTURE (IDL SAVE FILE) AND RUN TOMO_MAKE_GEOM, SAVING TO GEOM.PKL IN DIRECTORY ABOVE
if os.path.exists(savenamegeom) and not overwrite_geom:
    print("Reading ",savenamegeom)
    afile=open(savenamegeom,"rb")
    geom=pickle.load(afile)
    afile.close()
else:
    geom=tomo.tomo_make_geom(d)
    afile=open(savenamegeom,"wb")
    pickle.dump(geom,afile)
    afile.close()

if os.path.exists(savenamesphrecon) and not overwrite_sphrecon:
    print("Reading ",savenamesphrecon)
    afile=open(savenamesphrecon,"rb")
    sphrecon=pickle.load(afile)
    afile.close()
    print("Reading ",savenamesphdata)
    afile=open(savenamesphdata,"rb")
    sphdata=pickle.load(afile)
    afile.close()
else:
    sphrecon,sphdata=tomo.tomo_prep_sph(geom)
    afile=open(savenamesphrecon,"wb")
    pickle.dump(sphrecon,afile)
    afile.close()
    afile=open(savenamesphdata,"wb")
    pickle.dump(sphdata,afile)
    afile.close()

if os.path.exists(savefit) and not overwrite_fit:
    print("Reading ",savefit)
    afile=open(savefit,"rb")
    t=pickle.load(afile)
    afile.close()
else:
    t=tomo.tikhonov_search_tomo(d,sphrecon,sphdata,geom)
    afile=open(savefit,"wb")
    pickle.dump(t,afile)
    afile.close()

if os.path.exists(savecoeff) and not overwrite_coeff:
    print("Reading ",savecoeff)
    afile=open(savecoeff,"rb")
    c=pickle.load(afile)
    afile.close()
else:
    c=tomo.tikhonov_opt_tomo(t)
    afile=open(savecoeff,"wb")
    pickle.dump(c,afile)
    afile.close()

sh=np.shape(sphdata["sph"])
npa=sh[0]
nt=sh[1]
nsph=sh[2]
bmod=np.sum(np.broadcast_to(c,(npa,nt,nsph))*sphdata["sph"],axis=2)
sh=np.shape(sphrecon["sph"])
nlon=sh[0]
nlat=sh[1]
nsph=sh[2]
dens=np.sum(np.broadcast_to(c,(nlon,nlat,nsph))*sphrecon["sph"],axis=2)

plt.subplot(3,1,1)
plt.title("Model pB")
plt.imshow(np.transpose(bmod),origin='lower')
plt.subplot(3,1,2)
plt.title("Observed pB")
plt.imshow(d["im"],origin='lower')
plt.subplot(3,1,3)
plt.title("Model density")
plt.imshow(np.transpose(dens),origin='lower')
plt.savefig('cortom.png')

plt.close()


