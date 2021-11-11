import numpy as np
import math
import general_huw as gen
import time_huw
import coordinates_huw as coord
import constants_huw as const
import astropy.time as Time
from sunpy.time import parse_time
from stereo_spice.coordinates import StereoSpice

def makeg(p,r,x={},y={},z={},nocube=False,bk=False,nopos=False):

    if isinstance(x,dict) == False:
        if nocube == False:
            nx=np.size(x)
            ny=np.size(y)
            nz=np.size(z)
            xx, yy, zz = np.meshgrid(x,y,z)
            r=np.sqrt(xx**2+yy**2+zz**2)
            if nopos == True:
                p=1
            else:
                p=np.sqrt(yy**2+zz**2)
        else:
            r=np.sqrt(x**2+y**2+z**2)
            if nopos == True:
                p=1
            else:
                p=np.sqrt(y**2+z**2)

    s=1/r
    s2=s**2
    c=np.sqrt(1-s2)
    ar=c*s2
    br=-(1-3*s2-(c**2)*((1+3*s2)/s)*np.log((1+s)/c))/8.
    u=0.5607 #limb_dark_c2('orange')
    tsc=7.95e-26 #correct
    # tsc = 1.24878D-25; constant is this unless
    con=math.pi*tsc*0.5/(1-(u/3.))

    if bk == False:
        g=(p**2)*((1-u)*ar+u*br)/(r**2)
    else:
        g0=(p**2)*((1-u)*ar+u*br)/(r**2)
        cr=(4/3.)-c-((c**3)/3.)
        dr=(5+s2-(c**2)*((5-s2)/s)*np.log((1+s)/c))/8.
        g1=(1-u)*cr+u*dr
        g=2*g1-g0

    g=con*g

    return g

def tomo_make_geom(d,rmain,nt,npa,nx,dates={},spacecraft={}):

    if isinstance(spacecraft,dict) == True:
        spacecraft=d["system"]

    if isinstance(dates,dict) == True:
        dates=d["dates"]

    n=np.size(dates)
    if nt == n:
        tai=time_huw.anytim2tai(dates)
    else: 
        tai=gen.congrid(time_huw.anytim2tai(dates),nt)

    dates=time_huw.anytim2cal(tai,form=11,tai=True,msec=False)

    pa=coord.make_coordinates(npa,[0,360],minus_one=True)*np.deg2rad(1)

    rsun_cm=const.phys_constants(rsun=True)

    print('Calculating geometry...')
    spice = StereoSpice()
    obstime=parse_time(dates)
    distlonlat=spice.get_lonlat(obstime,spacecraft,'CARRINGTON')
    dist=distlonlat[0:,0]/(rsun_cm*1.e-5)

    xobsfact=rmain*2.5
    xobs=coord.make_coordinates(nx,[-1,1])*xobsfact

    lonsc,dist,xxobs = np.meshgrid(pa,dist,xobs,indexing="ij")
    latsc=np.arctan2(rmain,dist)
    dx=(xobs[1]-xobs[0])*rsun_cm #is the same for all LOS
    rsc=xxobs+dist; # distance from spacecraft to each LOS point;sqrt(dist^2-rmain^2)

    # convert to Heliocentric-Cartesian (section 3.1 Thompson 2006 https://www.aanda.org/articles/aa/pdf/2006/14/aa4262-05.pdf)
    y,x,z = coord.spherical2cartesian(rsc,lonsc,latsc)#y,x,z order because position angle (lonsc) measured from north
    z=dist-z
    x=-x # because position angle (lonsc) measured CCW from north
    x0=np.copy(x)
    y0=np.copy(y)
    z0=np.copy(z)

    # convert to Stonyhurst (section 7 of Thompson 2006), 
    # and conversion to Carrington (longitude only for this, see Thompson eq. 3)
    # Wish I could just use the StereoSpice conversion routines directly (like I do in IDL), 
    # but convert_coord does not accept "HGRTN" with the observer as STEREO A, it expects observer as Sun.
    # May contact Luke re. this 2021/11.
    carrearth=spice.get_lonlat(obstime,"Earth",'CARRINGTON')
    l0=carrearth[:,1]
    loncarrobs=distlonlat[:,1]
    lonstobs=np.deg2rad(loncarrobs-l0)
    b0=np.deg2rad(distlonlat[:,2])
    r=np.sqrt(x**2+y**2+z**2)
    lat=np.zeros((npa,nt,nx))
    lon=np.zeros((npa,nt,nx))
    for i in range(nt):
        lat[:,i,:]=np.arcsin((y[:,i,:]*np.cos(b0[i])+z[:,i,:]*np.sin(b0[i]))/r[:,i,:])
        lonsh=lonstobs[i]+np.arctan2(x[:,i,:],z[:,i,:]*np.cos(b0[i])-y[:,i,:]*np.sin(b0[i]))
        lon[:,i,:]=gen.wrap_n(np.rad2deg(lonsh)+l0[i],360)

    colat=np.deg2rad(90)-lat
    lon=np.deg2rad(lon)
    # convert to cartesian. x,y,z are now cartesian Carrington coordinates
    x,y,z=coord.spherical2cartesian(r,lon,colat)

    g=makeg(rmain,r,bk=(spacecraft=="soho"))
    # ***variables all agree with IDL up to this point***

    rdropoffpwr=2.2
    rdropoff=(rmain**rdropoffpwr)/(r**rdropoffpwr)
    geomult=g*dx*rdropoff
    tot_g=np.sum(geomult,axis=2)

    # calculate correction factor
    xc=coord.make_coordinates(1001,[-1,1])*rmain*3
    rc=np.sqrt(xc**2+rmain**2)
    gc=makeg(rmain,rc,bk=(spacecraft=='soho'))
    rdropoffc=(rmain**rdropoffpwr)/(rc**rdropoffpwr)
    dxc=(xc[1]-xc[0])*rsun_cm
    geomultc=gc*dxc*rdropoffc
    tot_gc=np.sum(geomultc)
    corr_fact=tot_gc/np.mean(tot_g)
    tot_g=tot_g*corr_fact

    geom = {
        "dates":dates,
        "pa":pa,
        "rmain":rmain,
        "x":x,
        "y":y,
        "z":z,
        "g":g,
        "dx":dx,
        "tot_g":tot_g,
        "r":r,
        "lon":lon,
        "lat":colat,
        "geomult":geomult,
        "corr_fact":corr_fact
    }
    
    return geom

