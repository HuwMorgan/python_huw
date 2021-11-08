import numpy as np
import matplotlib.pyplot as plt
import filters_huw as filt

a=np.random.normal(0,1,100)
a[10]=10
a[50]=11
a[76]=-10
plt.plot(a,'b+',label='Unfiltered')
b=filt.sliding_window_1d(a,9)
md=filt.sliding_window_1d(a,9,type="median")
mn=filt.sliding_window_1d(a,9,type="min")
mx=filt.sliding_window_1d(a,9,type="max")
#mx=filt.max_1d(a,9)
#mn=filt.min_1d(a,9)
plt.plot(b,'r',label='Mean')
plt.plot(md,'r--',label='Median')
plt.plot(mx,'y',label='Max')
plt.plot(mn,'y--',label='Min')
c=filt.point_filter_1d(a,wdth=5,median=True,sig=2.5)
plt.plot(c,'g',label='Point filter')
plt.legend()
plt.show()