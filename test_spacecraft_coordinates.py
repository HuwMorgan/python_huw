import astropy.time as Time
from stereo_spice.coordinates import StereoSpice

spice = StereoSpice()
files = spice.__get_kernal_files__()

obstime=Time("2017-08-01T12:00:00")
sta=spice.get_lonlat(obstime,'sta','CARRINGTON',precess=False)

a=sta