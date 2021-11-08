import os.path as path
import numpy as np
import time_huw 
import files_huw as files

file=['/Users/hum2/data/stereo/secchi/pb/a/cor2/2017/07/07/cor2_a_pb_20170707060815.dat', \
        '/Users/hum2/data/stereo/secchi/pb/a/cor2/2017/07/07/cor2_a_pb_20170707060815.dat']

date=files.filename2date_huw(file,tai=True)
print(file)
print(time_huw.anytim2cal(date,form=11))