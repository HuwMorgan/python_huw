import numpy as np

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