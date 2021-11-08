import os
import numpy as np
import matplotlib.pyplot as plt
import functions_huw as func
import coordinates_huw as coord

print('Testing spherical2cartesian and cartesian2spherical')
r=coord.make_coordinates(10,[2.2,8])
the=coord.make_coordinates(10,[10,350])
phi=coord.make_coordinates(10,[10,170])
x,y,z=coord.spherical2cartesian(r,the,phi,degrees=True,carrington=False)
r2,the2,phi2=coord.cartesian2spherical(x,y,z,degrees=True,wrap=True)
print('Ratio of output R to input R:', np.nanmin(r2/r),np.nanmax(r2/r))
print('Ratio of output Theta to input Theta:',np.nanmin(the2/the),np.nanmax(the2/the))
print('Ratio of output Phi to input Phi:',np.nanmin(phi2/phi),np.nanmax(phi2/phi))

nx=100
ny=70
xcen=nx/3
ycen=ny/2

g=func.gauss_2d(nx,ny,xcen=xcen,ycen=ycen,xsig=10,ysig=10)
g=np.where(g > 0.9,0,g)
x=np.arange(0,nx)-xcen
y=np.arange(0,ny)-ycen
# xx=np.ravel(np.resize(x,(nx,ny)),order='C')
# yy=np.ravel(np.swapaxes(np.resize(y,(ny,nx)),1,0),order='C')
xx=np.resize(x,(ny,nx))
yy=np.swapaxes(np.resize(y,(nx,ny)),1,0)
ang=np.arctan2(-xx/yy)
a=np.sin(ang*10)+np.sin(ang*5)+np.cos(ang*3)+np.cos(ang*2)
a=g*a

parange=[0,360]
rra=[5,np.sqrt(xcen**2+ycen**2)]
npa=300
nr=100
pix_size=1
b,pam,rm=coord.cartesian2polar(a,parange,rra,npa,nr,xcen,ycen,pix_size)

plt.subplot(211)
xyextent=[np.min(x),np.max(x),np.min(y),np.max(y)]
plt.imshow(a, cmap=plt.cm.BuPu_r,origin='lower',extent=xyextent,interpolation="none")
plt.subplot(212)
prextent=[parange[0],parange[1],rra[0],rra[1]]
plt.imshow(b, cmap=plt.cm.BuPu_r,origin='lower')

plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes([0.85, 0.1, 0.075, 0.8])
plt.colorbar(cax=cax)
filename='test_cartesian2polar.png'
plt.savefig(filename, bbox_inches='tight')
os.system('open '+filename)
# plt.show(block=False)
# plt.pause(5)
plt.close()
