import numpy as np
import math

# returns 'g', or geometrical weighting factor for Thomson scattering
# see Quemerais & Lamy 2002, A&A
# Calling: g=makeg(p,r)
# input p is plane-of-sky distance (heliocentric) of observation (e.g. 3 solar radii)
# input r is true heliocentric distance of point in corona (e.g. 5 solar radii)
# default is to return 'g' for polarized brightness (using orange filter of LASCO C2).
# if keyword set for Bk, e.g.  g=makeg(p,r,/bk), then the 'g' factor for K-coronal (total) brightness is returned. 

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
