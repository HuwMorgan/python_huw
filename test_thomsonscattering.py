import os
import numpy as np
import matplotlib.pyplot as plt
import functions_huw as func
import coordinates_huw as coord
import thomsonscattering_huw as thom

print('Testing makeg')

r=coord.make_coordinates(10,[2.2,8])
the=np.deg2rad(45)
p=r*np.sin(the)

g=thom.makeg(p,r)
g2=thom.makeg(p,r,bk=True)

plt.plot(g,'r',label='g for pB')
plt.plot(g2,'g',label='g for Bk')
plt.legend()

filename='test_makeg.png'
plt.savefig(filename, bbox_inches='tight')
os.system('open '+filename)
# plt.show(block=False)
# plt.pause(5)
plt.close()
