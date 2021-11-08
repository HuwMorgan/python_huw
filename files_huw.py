import os.path as path
import numpy as np
import time_huw 

def filename2date_huw(f00,tai=False,day=False):

    nf0=np.size(f00)

    if nf0 == 1 and not f00[0]:
        return f00[0]

    if isinstance(f00,list)==False:
        f0=[f00]
    else:
        f0=f00

    date=[]

    for if0 in range(nf0):
        
        f=path.basename(f0[if0])
        f=path.splitext(f)
        f=f[0]
  
        found_date=True

        fs0=f.split('_')
        fs=[fs0[i].split('.') for i in range(np.size(fs0))]
        nst=np.size(fs)
        if nst <= 1:
          print('files_huw.filename2date_huw: No dates found in filename (0) ',f0[if0])
          date.append('')
          found_date=False
        
        if found_date:

            isdate=np.empty(nst,dtype=np.int8)
            for i in range(nst):
                fnow=fs[i][0]
                for j in range(len(fnow)):
                    if fnow[j].isnumeric():
                        count=1
                    else:
                        count=0
                    isdate[i]=isdate[i]+count
      
        #   identify possible dates
            ind=np.where(isdate > 6)
            cnt=np.size(ind)
            if cnt==0:
                print('files_huw.filename2date_huw: No dates found in filename (1) ',f0[if0])
                found_date=False
                date.append("")

          
            if found_date:
                d=''
                ind=ind[0]
                for i in ind:
                    d=d+fs[i][0]
                date.append(time_huw.yyyymmdd2cal(d)[0])


    if day:
        date=time_huw.anytim2cal(date,form=11,date_only=True)
    if tai:
        date=time_huw.anytim2tai(date)

    if nf0==1 and isinstance(date,list)==True:
        date=date[0]


    return date