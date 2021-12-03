from sunpy.time import parse_time
from astropy.time import Time
import time_huw as thuw
import numpy as np

file='/Users'

date=['2011/10/01 11:57','2021/10/01 11:23']
print(date)
print()

print("Testing timegrid")
startdate="2011/01/01 00:00"
enddate="2011/01/10 00:00"
t=thuw.timegrid(startdate,enddate,delta=1,days=True)

print("Testing CAL:")
date2=thuw.anytim2cal(date,form=11)
print(date2)
print("")


print("Testing TAI:")
print(date)
tai=thuw.anytim2tai(date)
print(tai)
print("")

print("Testing form=8")
date11=thuw.anytim2cal(tai,form=8,date_only=False)
print(date11)
print("")

print("Testing yyyymmdd2cal:")
yyyymmdd=thuw.anytim2cal(date,form=8)
date3=thuw.yyyymmdd2cal(yyyymmdd)
print(yyyymmdd)
print(date3)
print("")