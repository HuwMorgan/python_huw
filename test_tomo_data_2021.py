import os
import numpy as np
import matplotlib.pyplot as plt
from tomo_data_2021 import tomo_data_2021


rmain=4.
middate='2017/07/07 12:00'

os.environ['PROCESSED_DATA'] = '/Users/hum2/data/sweep'
dir="/Users/hum2/data/sweep/"

d=tomo_data_2021(rmain, middate, overwrite=True)

plt.subplot(3,1,1)
plt.title("Model pB")
plt.imshow(np.transpose(d["bmod"]),origin='lower')
plt.subplot(3,1,2)
plt.title("Observed pB")
plt.imshow(np.transpose(d["bobs"]),origin='lower')
plt.subplot(3,1,3)
plt.title("Model density")
plt.imshow(np.transpose(d["dens"]),origin='lower')
plt.savefig(dir+'cortom.png')

plt.close()