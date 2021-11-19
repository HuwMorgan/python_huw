import inspect
import numpy as np
import astropy.wcs as wcs
from rebin import rebin
from scipy.interpolate import griddata
from scipy.interpolate import Rbf
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import RegularGridInterpolator
from .cor2_distortion import cor2_distortion
from ..display.movie.scc_sun_center import scc_sun_center

def bilinear_interp(XX, YY, ZZ, Xint, Yint):
    """ Bilinear interpolation a point in an unstructured spatial dataset.
        
        Usage
        -----
            >>> Zint = bilinear_interp(XX, YY, ZZ, Xint, Yint)
        
        Parameters
        ----------
            `XX`: array (float)
                X-cordinates of the known points
            `YY`: array (float)
                Y-cordinates of the known points
            `ZZ`: array (float)
                values of the known points (`XX`,`YY`,`ZZ` must have the same size)
            `Xint`: float, or vector/matrix of floats
                X-cordinate(s) of the point(s) to be interpolated
            `Yint`: float, or vector/matrix of floats
                Y-cordinate(s) of the point(s) to be interpolated
                `Xint` and `Yint` don't need to be structured 
                but must be the same size.
        
        Returns
        -------
            : float, or vector/matrix of floats
                the values to be interpolated
    """ 
    assert Xint.shape == Yint.shape
    rsize, csize = Xint.shape
    Zint = np.zeros((rsize,csize))

    for i in range(rsize):
        for j in range(csize):
            xpos, ypos = Xint[i,j], Yint[i,j]
            ip, jp = 0, 0
            # while XX[ip, jp] < xpos: # careful
            #       ip += 1 

            Zint[i,j] = Xint[i,j] + Yint[i,j]
    
    return Zint

def cor2_warp_huw(image, hdr, INFO=[]):

    if 'OBSRVTRY' in hdr:
        # determine key variables
        # w = wcs.WCS(hdr)
        # wcs = fitshead2wcs(hdr, system='A')
        sc = hdr['OBSRVTRY']
        sc = (''.join(sc.split()))[-1]
        # sc = strcompress(sc, remove_all=True)[-1]
        # Find corresponding control points x0 and y0 from x and y: 
        # Warp is determined as a function of sun's center.
        
        # Finds sun center in pixels.  Note that scc_sun_center returns the
        # sun center based on information in the image header# therefore, this
        # value has been scaled based on the image size.  In other words, for
        # beacom data, suncen=[128,128] (approximately)
        suncen = scc_sun_center(hdr)
        xc = suncen['xcen']
        yc = suncen['ycen']

        # Calculate the amount of binning that has been applied to the input
        # image# a similar amount of binning must shortly be applied to the
        # control points.
        sumxy = 2**(hdr['summed'] - 1) #;if image is 2048x2048, sumxy is 1
                                        #;if image is 1024x1024, sumxy is 2

    else:
        # get values from mvi header
        pos = hdr['filename'].find('.ft')
        sc = hdr['filename'][pos-1:pos]
        xc = hdr['xcen']
        yc = hdr['ycen']
        sumxy = hdr['cdelt1']/14.7    #14.7 arcsec/pixel for full res COR2 A and B


  
    #;query size of image
    sz=np.shape(image)
    nx=sz[0]
    ny=sz[1]
    
  
    #;create 2D arrays containing x indices [0,1,2,3,..,nx-1] in rows and y indices [0,1,2,3,..,ny-1] in columns
    #;correct for possible image rebinning with sumxy
    x=np.arange(nx)/sumxy
    y=np.arange(ny)/sumxy
    xout , yout = np.meshgrid(x,y)
    rout = np.sqrt((xout-sumxy*xc)**2+(yout-sumxy*yc)**2)  # calculate distance of each pixel from image center

#   ;the image warp is based on distance r from image center
#   ;cor2_distortion takes the input r for each pixel and specifies the new distance rout
#   ;the pixel value should be moved to (warp image from r to rout)
    
    rreg=np.linspace(0,np.amax(rout)*1.01,num=1001)
    routreg,cf,histinfo= cor2_distortion(rreg, sc)#[0] / sumxy
    rin = np.interp(rout,routreg,rreg) # rin agrees with IDL

    theta=np.arctan2((yout-yc*sumxy),(xout-xc*sumxy))# as above, agrees with IDL
    xin= np.swapaxes((rin*np.cos(theta)+xc*sumxy)/sumxy,1,0)
    yin= np.swapaxes((rin*np.sin(theta)+yc*sumxy)/sumxy,1,0)

    image=np.where(np.isfinite(image),image,0.)
    interp_spline = RectBivariateSpline(x, y, image)
    outimage = interp_spline.ev(xin, yin) # and this agrees with IDL!
  
    id = '$Id: cor2_warp_huw.py,v 1.0 2021/11 hmorgan@aber.ac.uk $'
    histinfo = [INFO, id[1:-2]]
  
    return outimage , histinfo
