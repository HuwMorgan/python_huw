import numpy as np
import math

def congrid(a,n_out):

    # equivalent of IDL's congrid with keywords /interp and /minus_one

    n_in=np.size(a)
    if n_in == n_out:
        return a
    
    xout=np.linspace(0,n_in-1,num=n_out)
    xin=np.arange(0,n_in)
    b = np.interp(xout,xin,a)

    return b

def minmax(a):
    min=np.nanmin(a)
    max=np.nanmax(a)
    mnmx=[min,max]
    return mnmx

def pad_2d(a,n=1):

    sh=np.shape(a)
    nx=sh[0]
    ny=sh[1]
    a2=np.zeros(nx+n*2,ny+n*2)
    a2[1,1]=a

    return a2

def rebin(a,n_out):

    print("rebin: this function is not yet written properly, use congrid for interpolation of 1D vectors")

    n_in=np.size(a)
    if n_in == n_out:
        return a

    rat=n_in/n_out
    print("Ratio = ",rat)
    if rat < 1:

        rat=n_out/n_in
        step= n_out//n_in
        
        if step != rat:
            print("rebin: n_out should be integer factor of n_in!")
            return
        
        xout=np.linspace(0,n_in-1,num=n_out)
        xin=np.arange(0,n_in)
        b = np.interp(xout,xin,a)

    else:
        step= n_in//n_out
        
        if step != rat:
            print("rebin: n_out should be integer factor of n_in!")
            return

        b = np.empty(n_out) 
        for iout,iin in enumerate(np.arange(0,n_in,step)):
            b[iout] = np.mean(a[iin:iin+step])

    return b

def sizearr(a):
    ndim=np.ndim(a)
    sh=np.shape(a)
    if ndim == 1:
        return sh[0]
    if ndim == 2:
        return sh[0],sh[1]
    if ndim == 3:
        return sh[0],sh[1],sh[2]
    if ndim == 4:
        return sh[0],sh[1],sh[2],sh[3]
    if ndim == 5:
        return sh[0],sh[1],sh[2],sh[3],sh[4]
        
def wrap_n(a,wrap):
    b=((a % wrap)+wrap) % wrap
    return b
