import numpy as np
from scipy.ndimage.filters import uniform_filter1d
from scipy import signal
import functions_huw as func

#sliding window mean/median/max/min smoothing for 1D vectors, default is mean smoothing.
#x is input vector
#wdth is width of sliding window (default=9)
#type can be "mean" (default), "median", "max", or "min"
#Note that towards edges of input vector, the sliding window is truncated
#Huw Morgan at Aberystwyth University, 2021/10, hmorgan@aber.ac.uk
def sliding_window_1d(x,wdth=9,type="mean"):

    lenx=len(x)
    if lenx <= 1:
        print("x has <2 elements, returning x (sliding_window_1d)")
        return x

    n2=wdth//2
    xm=np.empty(lenx)
    for i in range(lenx):
        
        ileft=i-n2
        
        if ileft < 0:
            ileft=0
        
        iright=i+n2
        if iright > (lenx-1):
            iright=lenx-1
        
        #want to use match (case) statement here, but only available on python 3.10.
        if type.lower()=="mean":
            xm[i]=np.mean(x[ileft:iright])
        if type.lower()=="median":
            xm[i]=np.median(x[ileft:iright])
        if type.lower()=="min":
            xm[i]=np.min(x[ileft:iright])
        if type.lower()=="max":
            xm[i]=np.max(x[ileft:iright])
            # case _:
            #     print("sliding_window_1d type not recognised, returning original vector")
            #     xm=x
            #     break
        #print(i,ileft,iright,lenx)

    return xm


#For 1D vectors, identifies outliers by comparing values to local standard deviation. Replaces value by new value closer to local mean. Iterates.
#c0 is input vector
#wdth is width of sliding window (default=9)
#nsig is significance value, or multiplier of standard deviation that sets threshold to identify outliers (default 1.6 ~80%)
#niter is maximum number of iterations (default=12)
#log (default False), if set to True, applies algorithm to log10 of input vector. Note, converts back to non-log after algorithm
#median (default False), if set to True, then instead of local mean and standard deviation, uses local median and local median absolute deviation.
#silent (default False), if set to True, then suppresses print at each iteration
#difffactor is factor by which outlying values are changed closer to the local mean (default 0.5). 
#So new value of outlier will be set at halfway between local mean and outlying value.
#Huw Morgan at Aberystwyth University, 2021/10, hmorgan@aber.ac.uk
def point_filter_1d(c0,wdth=9,sig=1.6,niter=12, 
    log=False, median=False, silent=False, difffactor=0.5):

    if len(c0) < wdth:
        print("Point_filter_1d: length of signal shorter than smoothing window")
        return c0

    if log:
        c=np.log10(c0) 
    else:
        c=c0
    
    cnt=1
    iter=0
    if silent == False:
        print('Iteration','N outliers')

    while (cnt > 0) and (iter < niter):

        if median:
            cs=sliding_window_1d(c,wdth=wdth)
        else:
            cs=uniform_filter1d(c,size=wdth,mode='reflect')
        
        df=c-cs
        if median:
            st=sliding_window_1d(abs(df),wdth=wdth)
        else:
            st=np.sqrt(uniform_filter1d(df**2,size=wdth,mode='reflect'))

        cnt=np.count_nonzero(abs(df)>(st*sig))
        if cnt > 0:
            ind=np.nonzero(abs(df)>(st*sig))
            c[ind]=difffactor*df[ind]+cs[ind]
        
        if silent == False:
            print(iter,'      ',cnt)
        
        iter=iter+1
    
    if log:
        c=10**c
    
    return c


#For 2D arrays, identifies outliers by comparing values to local standard deviation. Replaces value by new value closer to local mean. Iterates.
#c0 is input array
#width is width of sliding window (default=3)
#nsig is significance value, or multiplier of standard deviation that sets threshold to identify outliers (default 4)
#niter is maximum number of iterations (default=12)
#silent (default False), if set to True, then suppresses print at each iteration
#difffactor is factor by which outlying values are changed closer to the local mean (default 0.5). 
#So new value of outlier will be set at halfway between local mean and outlying value.
#Huw Morgan at Aberystwyth University, 2021/10, hmorgan@aber.ac.uk
def point_filter_2d(im0,width=3,nsig=2.6,niter=8, 
    silent=False, difffactor=0.5):

    if np.size(width) == 1:
        w=[width,width]
    else:
        w=width

    kxy=func.gauss_2d(w[0]*10,w[1]*10,xsig=w[0],ysig=w[1])
    ind=np.nonzero(kxy>0.2)
    ind0=np.amin(ind,axis=1)
    ind1=np.amax(ind,axis=1)
    kxy=kxy[ind0[0]:ind1[0],ind0[1]:ind1[1]]
    kxy=kxy/np.sum(kxy)
    nkxy=np.shape(kxy)

    nxy=np.shape(im0)
    nx=nxy[0]
    ny=nxy[1]
    if (nkxy[0] > nx) or (nkxy[1] > ny):
        print("filters_huw.point_filter_2d: one or more image dimensions too small, returning")
        return im0

    im=np.copy(im0)
    iter=0
    cnt=1.e6

    if silent == False:
        print('Iteration | ','N outliers')

    while (cnt > 10) and (iter < niter):

        av=signal.convolve2d(im,kxy,boundary='symm',mode='same')
        
        df=im-av
        st=np.sqrt(signal.convolve2d(df**2,kxy,boundary='symm',mode='same'))

        cnt=np.count_nonzero(abs(df)>(st*nsig))
        if cnt > 0:
            ind=np.nonzero(abs(df)>(st*nsig))
            im[ind]=av[ind]+difffactor*df[ind]
        
        if silent == False:
            print(iter,'        | ',cnt)
        
        iter=iter+1
    
    return im
