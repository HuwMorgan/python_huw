import numpy as np
import scipy.special as spspec
import math

def gauss_2d(nx,ny,xcen={},ycen={},xsig=1,ysig=1,norm=False):

    # Initializing values of x-axis and y-axis
    if  isinstance(xcen,dict) == True:
        xcen=nx2=nx//2
    if  isinstance(ycen,dict) == True:
        ycen=ny2=ny//2
    xx, yy = np.meshgrid(np.arange(0,nx)-xcen, np.arange(0,ny)-ycen)
    
    # Calculating Gaussian array
    gauss = np.exp(-( (xx**2) / ( 2.0 * xsig**2 ) )-( (yy**2) / ( 2.0 * ysig**2 ) )) 

    if norm == True:
        gauss=gauss/np.sum(gauss)
    
    return gauss

def make_spher_harm(nl,lon_in,lat,userlm=None):

    pi=math.pi
    lon=np.where(lon_in < 0,lon_in+2*pi,lon_in)

    if userlm == None:
        nsph=nl*(nl+1)+nl+1
        lstend=[0,nl]
    else:
        nsph=1
        lstend=[userlm[0],userlm[0]]
    
    # ;note factor 'f' below (compare IDL help for spher_harm and e.g. equations B.3 of appendix B
    # ;to do with use of real-valued spherical harmonics at all m. At m=0, imaginary part is zero everywhere so
    # ;sqrt(2) factor not valid/needed.

    ndim=np.ndim(lon)
    sz=np.shape(lon)
    if ndim == 0: 
        sph=np.zeros(nsph)
    else:
        sph=np.zeros(np.append(sz,nsph))
    
    isph=0
    for l in range(lstend[0],lstend[1]+1):
        mstend=[-l,l] if userlm == None else [userlm[1],userlm[1]]
        for m in range(mstend[0],mstend[1]+1):
            f=1 if m == 0 else np.sqrt(2)
            sphnow=spspec.sph_harm(m,l,lon,lat)
            sphnow=np.imag(sphnow) if m<0 else np.real(sphnow)
            sphnow=sphnow*f
            if ndim == 0:sph[isph]=sphnow
            if ndim == 1:sph[:,isph]=sphnow
            if ndim == 2:sph[:,:,isph]=sphnow
            if ndim == 3:sph[:,:,:,isph]=sphnow # yes there are more elegant ways of doing this...!
            isph+=1
    
    if nsph==1: sph=np.squeeze(sph)

    return sph