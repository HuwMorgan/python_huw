import numpy as np
import math
import general_huw as gen
import time_huw as time

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

def tomo_make_geom(d,rmain,nt,npa,nx,dates={},spacecraft={}):

    if isinstance(spacecraft,dict) == True:
        spacecraft=d["system"]

    if isinstance(dates,dict) == True:
        dates=d["dates"]

    n=np.size(dates)
    if nt == n:
        tai=time.anytim2tai_huw(dates)
    else: 
        tai=gen.congrid_huw(time.anytim2tai_huw(dates),nt)

    dates=time.anytim2cal_huw(tai,form=11)

dates=anytim2cal(tai,form=11)
pa=make_coordinates(npa,[0,360],/minus)*!dtor
tt=(rebin(reform(tai,1,nt),npa,nt))
ppa=(rebin(pa,npa,nt))
rsun_cm=phys_constants(/rsun)

print,'Calculating geometry...'

dist=get_stereo_lonlat(dates,spacecraft,system='Carr',/deg)
dist=reform(dist[0,*])/(rsun_cm*1.e-5)
dist=rebin(reform(dist,1,nt),npa,nt,nx)

xobs=make_coordinates(nx,[-1,1])
;xobs=xobs*rmain*2

xobs=rebin(reform(xobs,1,1,nx),npa,nt,nx)
xobsfact=rmain*2.5
xobs=temporary(xobs)*xobsfact;rebin(xobsfact,npa,nt,nx)
dx=rebin((xobs[*,*,1]-xobs[*,*,0])*rsun_cm,npa,nt,nx)
rsoho=xobs+dist;distance from SOHO to each LOS point;sqrt(dist^2-rmain^2)
latsoho=atan(rmain/dist);asin(rmain/rsoho) HUW 2014/10/29
lonsoho=rebin(ppa,npa,nt,nx)

;convert to cartesian and to HGRTN
spherical2cartesian,rsoho,lonsoho,latsoho,z,y,x

x=dist-temporary(x);heliocentric
y=-temporary(y)
r=sqrt(x^2+y^2+z^2)

g=fltarr(npa,nt,nx)

for i=0,nt-1 do begin
  if i mod 50 eq 0 then print,i,' out of ',nt-1
  coord=rotate([[(x[*,i,*])[*]],[(y[*,i,*])[*]],[(z[*,i,*])[*]]],4)
  coord2=coord
  convert_stereo_coord,dates[i],coord,'HGRTN','Carr',spacecraft=spacecraft
  x[*,i,*]=coord[0,*]
  y[*,i,*]=coord[1,*]
  z[*,i,*]=coord[2,*]
  g[*,i,*]=makeg(rmain,r[*,i,*],bk=spacecraft eq 'soho')
endfor

cartesian2spherical,x,y,z,r,lon,lat

rdropoffpwr=2.2
rdropoff=(rmain^rdropoffpwr)/(r^rdropoffpwr)
geomult=g*dx*rdropoff
tot_g=total(geomult,3)

;calculate correction factor
xc=make_coordinates(1001,[-1,1])*rmain*3
rc=sqrt(xc^2+rmain^2)
gc=makeg(rmain,rc,bk=spacecraft eq 'soho')
rdropoffc=(rmain^rdropoffpwr)/(rc^rdropoffpwr)
dxc=(xc[1]-xc[0])*rsun_cm
geomultc=gc*dxc*rdropoffc
tot_gc=total(geomultc)
corr_fact=tot_gc/mean(tot_g)
tot_g=tot_g*corr_fact

geom={dates:dates,pa:pa,rmain:rmain,x:x,y:y,z:z,g:g, $
  dx:dx,tot_g:tot_g,r:r,lon:lon,lat:lat, $
  geomult:geomult,corr_fact:corr_fact}
  
return,geom

end
 
