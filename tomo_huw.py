import numpy as np
import math
import general_huw as gen
import functions_huw as func
import filters_huw as filt
import time_huw
import coordinates_huw as coord
import constants_huw as const
import astropy.time as Time
from sunpy.time import parse_time
from stereo_spice.coordinates import StereoSpice
from scipy.interpolate import RectBivariateSpline

import matplotlib.pyplot as plt


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

def tikhonov_opt_tomo(d):

    sh=np.shape(d["c"])
    nshp=sh[0]
    nl=sh[1]
    nd=sh[2]

    g=d["rho"]

    thresh=np.percentile(g,15)
    mask=g < thresh

    mask=filt.region_size_filter(mask,npix=5)

    nfine=512
    l=np.linspace(np.log10(np.nanmin(d["l"])),np.log10(np.nanmax(d["l"])),num=nfine)
    l=10**l
    dens=np.linspace(np.nanmin(d["dens"]),np.nanmax(d["dens"]),num=nfine)

    ix=np.interp(l,d["l"],np.arange(d["nl"]))
    iy=np.interp(dens,d["dens"],np.arange(d["ndens"]))
    ixx,iyy=np.meshgrid(ix,iy)
    interp = RectBivariateSpline(np.arange(d["nl"]), np.arange(d["ndens"]), g,kx=1,ky=1)
    gfine = interp.ev(ixx, iyy) 
    interpmsk = RectBivariateSpline(np.arange(d["nl"]), np.arange(d["ndens"]), mask,kx=1,ky=1)
    maskfine=interpmsk.ev(ixx,iyy)
    indf=maskfine > 0
    ixfine=np.arange(nfine)
    iyfine=np.arange(nfine)
    ixxfine,iyyfine=np.meshgrid(ixfine,iyfine)
    ixf=ixxfine[indf]
    iyf=iyyfine[indf]
    gf=gfine[indf]
    ixfineopt=np.sum(ixf*gf)/np.sum(gf)
    iyfineopt=np.sum(iyf*gf)/np.sum(gf)

    r=np.sqrt(ixf**2+iyf**2)
    indopt=np.nanargmax(r)
    ixfineopt2=ixf[indopt]
    iyfineopt2=iyf[indopt]

    ixfineopt=(ixfineopt+ixfineopt2)/2
    iyfineopt=(iyfineopt+iyfineopt2)/2

    ixopt=np.interp(ixfineopt,ixfine,ix)
    iyopt=np.interp(iyfineopt,iyfine,iy)
    lopt=np.interp(ixfineopt,ixfine,l)
    dopt=np.interp(iyfineopt,iyfine,dens)

    inv=np.linalg.inv(d["xx"]+lopt*d["id"])
    c=np.sum(d["xy2"]*inv,1)
    dens=np.sum(np.broadcast_to(c,(d["nlon"],d["nlat"],d["nsph"]))*d["sph"],axis=2)
    dens=np.where(dens > dopt,dens,dopt)
    c=[np.sum(np.squeeze(d["sph"][:,:,isph])*dens*d["sinlat"])*d["dlon"]*d["dlat"] for isph in np.arange(0,d["nsph"])]

    return c


def tikhonov_search_tomo(d,sphrecon,sphdata,geom,nl=25,ndens=20):

    bobs=d["im"]
    relnoise=d["relnoise"]

    sh=np.shape(sphdata["sph"])
    nsphdata=sh[2]
    nt=sh[1]
    npa=sh[0]
    n=npa*nt

    a=np.zeros((n,nsphdata))
    for i in np.arange(nsphdata):
        a[:,i]=np.ravel(sphdata["sph"][:,:,i],order='F')
    mna=np.mean(np.abs(sphdata["sph"]))
    a = a/mna

    y=np.ravel(bobs)/mna
    noise=np.ravel(relnoise)*y

    index=~np.isnan(y) & ~np.isnan(noise)
    y=y[index]
    noise=noise[index]
    n=np.size(y)
    a=a[index,:]

    weights=1/noise
    weights=weights/np.mean(weights)
    a=a*np.broadcast_to(np.expand_dims(weights,1),(n,nsphdata))
    y=y*weights

    xx=np.matmul(np.transpose(a),a)
    xy=[np.sum(a[:,i]*y) for i in np.arange(nsphdata)]
    xxinv=np.linalg.inv(xx)
    xy2=np.broadcast_to(xy,(nsphdata,nsphdata))

    w=np.sum(np.abs(sphdata["lm"]),axis=1)
    w=w/np.mean(w)
    id=np.identity(nsphdata)*np.broadcast_to(w,(nsphdata,nsphdata))

    diagxx=np.diag(xx)
    lmn=np.nanmin(diagxx)*0.1
    lmx=np.nanmax(diagxx)*2.0
    l=10**np.linspace(np.log10(lmn),np.log10(lmx),num=nl)

    minb=np.percentile(bobs,2)
    indmin=bobs < minb
    bmin=bobs[indmin]
    geotemp=np.swapaxes(np.sum(geom["geomult"],axis=2),1,0)
    geotemp=geotemp[indmin]
    densmin=bmin/geotemp
    min_d0=np.mean(densmin)*0.2
    max_d0=np.mean(densmin)*2.0
    mindens=np.linspace(min_d0,max_d0,num=ndens)
    
    sph=sphrecon["sph"]
    sh=np.shape(sph)
    nsph=sh[2]
    nlat=sh[1]
    nlon=sh[0]
    lon=sphrecon["lon"]
    lat=sphrecon["lat"]
    sinlat=np.broadcast_to(np.sin(lat),(nlon,nlat))
    dlon=np.median(lon[1:]-lon[0:-1])
    dlat=np.median(lat[1:]-lat[0:-1])

    rho=np.zeros((nl,ndens))
    cmain=np.zeros((nsphdata,nl,ndens))
    
    for il in np.arange(0,nl):
        print(il,' out of ',nl-1)
    
        inv=np.linalg.inv(xx+l[il]*id)
        c=np.sum(xy2*inv,1)

        dens=np.sum(np.broadcast_to(c,(nlon,nlat,nsph))*sph,axis=2)
        
        for idens in np.arange(0,ndens):
            densnow=np.where(dens > mindens[idens],dens,mindens[idens])
            c=[np.sum(np.squeeze(sph[:,:,isph])*densnow*sinlat)*dlon*dlat for isph in np.arange(0,nsph)]
            yf=np.sum(a*np.broadcast_to(np.expand_dims(c,0),(n,nsphdata)),axis=1)            
            rho[il,idens]=np.mean(np.sqrt((yf-y)**2)/noise)
            cmain[:,il,idens]=c

    t={
        "nl":nl,
        "ndens":ndens,
        "l":l,
        "dens":mindens,
        "rho":rho,
        "c":cmain,
        "meanabsa":mna,
        "xx":xx,
        "id":id,
        "xy2":xy2,
        "sph":sph,
        "nlon":nlon,
        "nlat":nlat,
        "nsph":nsph,
        "sinlat":sinlat,
        "dlon":dlon,
        "dlat":dlat
    }

    return t

def tomo_make_geom(d,large=True):

    rmain=d["rmain"]
    spacecraft=d["system"]
    dates=d["dates"]

    nx=151 if large else 41
    npa=d["npa"]
    dates=d["dates"]
    nt=np.size(dates)
    
    tai=time_huw.anytim2tai(dates)
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
    lon=np.where(lon > 180,lon-360,lon)
    lon=np.where(lon < -180,lon+360,lon)
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

def tomo_prep_sph(geom,nlon=540,nlat=270,norder=22):

    pi=math.pi

    lonh,lath,llonh,llath,dlon,dlat=coord.make_lonlat(nlon,nlat)

    npa,nt,nx=gen.sizearr(geom["lon"])

    print('Constructing spherical harmonic basis...')

    nl=norder
    nsph=nl*(nl+1)+nl+1
    lm=np.zeros((nsph,2))
    sphdens=np.zeros((nlon,nlat,nsph))
    sphmain=np.zeros((npa,nt,nsph))
    cenx_sph=np.zeros((npa,nt,nsph))
    xcen=nx//2

    s=geom["lon"]
    s=s[:,:,xcen]

    df=np.roll(s,-1,axis=1)-s
    db=s-np.roll(s,1,axis=1)
    dlon=0.5*(db+df)
    dlon[:,0]=df[:,0]
    dlon[:,nt-1]=db[:,nt-1]
    dlon=np.where(dlon > pi/2,dlon-pi,dlon)
    dlon=np.where(dlon < -pi/2,dlon+pi,dlon)
    
    isph=0
    for l in range(0,nl+1):
        for m in range(-l,l+1):
            lm[isph,0]=l
            lm[isph,1]=m
            sphdens[:,:,isph]=func.make_spher_harm(nl,llonh,llath,userlm=[l,m])
            sphnow=func.make_spher_harm(nl,geom["lon"],geom["lat"],userlm=[l,m])
            sphmain[:,:,isph]=np.sum(sphnow*geom["geomult"],axis=2)#LOS integrated SH
            cenx_sph[:,:,isph]=func.make_spher_harm(nl,geom["lon"][:,:,xcen],geom["lat"][:,:,xcen],userlm=[l,m])*geom["geomult"][:,:,xcen]/dlon
            isph+=1
            if (isph % 10)==0: print(isph," out of ",nsph-1)

    sphrecon={
        "nsph":nsph,
        "nl":nl,
        "coeff":np.zeros(nsph),
        "lm":lm,
        "lon":lonh,
        "lat":lath,
        "sph":sphdens
    }

    sphdata={
        "nsph":nsph,
        "nl":nl,
        "coeff":np.zeros(nsph),
        "lm":lm,
        "sph":sphmain,
        "cenx_sph":cenx_sph
    }

    return sphrecon,sphdata
