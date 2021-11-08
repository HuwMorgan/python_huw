import numpy as np
import math
from scipy.interpolate import RegularGridInterpolator
import general_huw as gen


def cartesian2polar(im0,parange,rra,npa,nr,xcen,ycen,pix_size=1,roll=0,missing=float("nan")):

    sz=np.shape(im0)
    nx=sz[1]
    ny=sz[0]
    x,y,ht,pa = get_ht_pa_2d(nx,ny,xcen,ycen,pix_size=pix_size)#note do not supply roll keyword here, 
                                                                #roll accounted for later in interpolation block

    pam=np.linspace(parange[0],parange[1],npa)
    rm=np.linspace(rra[0],rra[1],nr)
    ppam,rrm=np.meshgrid(pam,rm)

    x0=-rrm*np.sin(np.deg2rad(ppam))
    y0=rrm*np.cos(np.deg2rad(ppam))
    x00=x0*np.cos(np.deg2rad(roll))+y0*np.sin(np.deg2rad(roll))
    y00=-x0*np.sin(np.deg2rad(roll))+y0*np.cos(np.deg2rad(roll))

    x00=np.ravel(x00)
    y00=np.ravel(y00)

    ix=np.interp(x00,x,np.arange(0,nx))
    iy=np.interp(y00,y,np.arange(0,ny))

    int=RegularGridInterpolator((y,x),im0,bounds_error=False,fill_value=float('nan'))

    imp=np.empty(npa*nr)
    for i in np.arange(0,npa*nr):
        imp[i]=int([y00[i],x00[i]])
    #then need to reshape imp to (npa,nr)

    imp=np.reshape(imp,(nr,npa))

    return imp,pam,rm


def cartesian2spherical(x,y,z,degrees=False,wrap=False,carrington=False):

    r=np.sqrt(x**2+y**2+z**2)

    if carrington == True:
        the=np.arctan2(-y,x)
    else:
        the=np.arctan2(y,x)
    
    phi=np.arccos(z/r)

    if degrees==True:
        the=np.rad2deg(the)
        phi=np.rad2deg(phi)

    if wrap == True:
        if degrees==True:
            n=360.
        else:
            n=2*math.pi
        the=gen.wrap_n(the,n)
    
    return r, the, phi


def get_ht_pa_2d(nx,ny,xcen,ycen,pix_size=1,roll=0):

    x=(np.arange(0,nx)-xcen)*pix_size
    y=(np.arange(0,ny)-ycen)*pix_size
    xx, yy = np.meshgrid(x,y)

    ht=np.sqrt(xx**2+yy**2)
    pa=np.rad2deg(np.arctan2(-xx,yy))+roll
    pa=np.where(pa > 360,pa-360,pa)
    pa=np.where(pa < 0,pa+360,pa)

    return x,y,ht,pa


def make_coordinates(n,range,minus_one=False,reverse=False):

    if n == 1:
        return range[0]

    if minus_one == True:
        endpoint=False
    else:   
        endpoint=True

    c=np.linspace(range[0],range[1],num=n,endpoint=endpoint)
    
    if reverse == True:
        c=np.flip(c)

    return c

def spherical2cartesian(r,the,phi,degrees=False,carrington=False,grid=False):

    if degrees==True:
        factor=np.deg2rad(1)
    else:
        factor=1
    
    nr=np.size(r)
    nthe=np.size(the)
    nphi=np.size(phi)
    if grid == False and (nthe != nphi) and (nthe > 1 and nphi > 1):
        print('Number of the and phi elements should be equal if grid keyword not set!')
        return

    if grid == True:
        pass
        # tthe=rebin(the,nthe,nphi)
        # pphi=rebin(reform(phi,1,nphi),nthe,nphi)
        # x=r*np.cos(tthe*factor)*np.sin(pphi*factor)
        # y=r*np.sin(tthe*factor)*np.sin(pphi*factor)
        # z=r*np.cos(pphi*factor)
    else:
        x=r*np.cos(the*factor)*np.sin(phi*factor)
        y=r*np.sin(the*factor)*np.sin(phi*factor)
        z=r*np.cos(phi*factor)

    if carrington == True:
        y=-1*y
        
    return x, y, z




