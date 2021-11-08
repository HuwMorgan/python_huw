from datetime import datetime
from sunpy.time import parse_time
# from sunpy.time import TimeTaiSeconds
from astropy.time import Time
import time_huw as thuw

file='/Users'

date=['2011/10/01 11:57','2021/10/01 11:23']
print(date)
print()

print("Testing CAL:")
date2=thuw.anytim2cal(date,form=11)
print(date2)
print("")


print("Testing TAI:")
print(date)
tai=thuw.anytim2tai(date)
print(tai)
print("")

print("Testing form=11")
date11=thuw.anytim2cal(tai,form=8,date_only=False)
print(date11)
print("")

print("Testing yyyymmdd2cal:")
yyyymmdd=thuw.anytim2cal(date,form=8)
date3=thuw.yyyymmdd2cal(yyyymmdd)
print(yyyymmdd)
print(date3)
print("")