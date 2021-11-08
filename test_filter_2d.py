import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import filters_huw as filt
import functions_huw as func

nx=100
ny=100
a=np.random.normal(0,1,size=(nx,ny))
b=filt.point_filter_2d(a,difffactor=0.)
a[10,90]=10
a[50,70]=-10
a[80,10]=12
print('Non filtered outlying pixel value =',a[10,90])
print('Point-filtered outlying pixel value =',b[10,90])
# f,axarr = plt.subplots(2)
# axarr[0].imshow(a)
# axarr[1].imshow(b)
# plt.show()
plt.subplot(211)
plt.imshow(a, cmap=plt.cm.BuPu_r)
plt.subplot(212)
plt.imshow(b, cmap=plt.cm.BuPu_r)

plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes([0.85, 0.1, 0.075, 0.8])
plt.colorbar(cax=cax)
plt.colorbar.clim(vmin, vmax) 
plt.show()