import numpy as np
import math

def congrid_huw(a,n_out):

    # equivalent of IDL's congrid with keywords /interp and /minus_one

    n_in=np.size(a)
    if n_in == n_out:
        return a
    
    xout=np.linspace(0,n_in-1,num=n_out)
    xin=np.arange(0,n_in)
    b = np.interp(xout,xin,a)

    return b

def rebin_huw(a,n_out):

    print("rebin_huw: this function is not yet written properly, use congrid for interpolation of 1D vectors")

    n_in=np.size(a)
    if n_in == n_out:
        return a

    rat=n_in/n_out
    print("Ratio = ",rat)
    if rat < 1:

        rat=n_out/n_in
        step= n_out//n_in
        
        if step != rat:
            print("rebin_huw: n_out should be integer factor of n_in!")
            return
        
        xout=np.linspace(0,n_in-1,num=n_out)
        xin=np.arange(0,n_in)
        b = np.interp(xout,xin,a)

    else:
        step= n_in//n_out
        
        if step != rat:
            print("rebin_huw: n_out should be integer factor of n_in!")
            return

        b = np.empty(n_out) 
        for iout,iin in enumerate(np.arange(0,n_in,step)):
            b[iout] = np.mean(a[iin:iin+step])

    return b

def wrap_n(a,wrap):
    b=((a % wrap)+wrap) % wrap
    return b