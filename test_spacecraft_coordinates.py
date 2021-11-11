import astropy.time as Time
from sunpy.time import parse_time
from stereo_spice.coordinates import StereoSpice

spice = StereoSpice()
obstime=parse_time("2017/08/01 12:00:00")
sta=spice.get_lonlat(obstime,'sta','CARRINGTON',precess=False)
print(sta)
