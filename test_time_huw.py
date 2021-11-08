from datetime import datetime
from sunpy.time import parse_time
# from sunpy.time import TimeTaiSeconds
from astropy.time import Time
import time_huw as thuw

date=['2011/10/01 11:57','2021/10/01 11:23']
date2=parse_time(date)

print("Testing TAI:")
print(date)
tai=thuw.anytim2tai_huw(date)
print(tai)
print("")


print("Testing form=11")
date11=thuw.anytim2cal_huw(tai,form=8,date_only=False)
print(date11)