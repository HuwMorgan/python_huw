import numpy as np
import matplotlib.pyplot as plt
import filters_huw as filt
import os
from scipy.io import readsav
import pickle

dir="/Users/hum2/data/cortomtests/pb_calibrate_secchi"
filepy= dir+"/cor2_a_pb_20170707000815.pkl"
fileidl=dir+"/cor2_a_pb_20170707000815.dat"

afile=open(filepy,"rb")
dpy=pickle.load(afile)
afile.close()

# FOLLOWING TO OPEN TOMO DATACUBE STRUCTURE (IDL SAVE FILE) AND RUN TOMO_MAKE_GEOM, SAVING TO GEOM.PKL IN DIRECTORY ABOVE
d={}
v=readsav(fileidl,verbose=False)
tag_name=list(v.keys())
values=list(v.values())
for i in range(np.size(tag_name)):
    print(tag_name[i])
    if tag_name[i]=="files":
        n=np.size(values[i])
        val = ["" for j in range(n)]
        for j in range(n):
            val[j]=values[i][j].decode("utf-8")
        valuenow=[val]
    elif tag_name[i]=="date" or tag_name[i]=="geometry":
        valuenow=values[i].decode("utf-8")
    else:
        valuenow=values[i]
    # if np.size(valuenow) > 1:
    #     valuenow=valuenow[0]
    d[tag_name[i]]=valuenow

print("Comparing pb_calibrate_secchi results")
print("Python pickle file | IDL save file")
print(dpy["date"],"|",d["date"])
print(dpy["para"],"|",d["para"])
print(dpy["rra"],"|",d["rra"])
print()

plt.plot(dpy["pb"][40,:],'b+',label='Python')
plt.plot(d["pb"][40,:],'r+',label='IDL')
plt.legend()
filename=dir+'/test_pb_calibrate_secchi.png'
plt.savefig(filename, bbox_inches='tight')
os.system('open '+filename)
# plt.show(block=False)
# plt.pause(5)
plt.close()