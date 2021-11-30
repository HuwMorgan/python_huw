from sunpy.time import parse_time
from astropy.time import Time
import numpy as np
import re

def anytim2tai(date_in):

    n=np.size(date_in)
    tai=np.empty(n)

    for i in range(n):
        if date_in[i] == "":
            tai[i]=float("nan")
            continue

        date=parse_time(date_in) if n==1 else parse_time(date_in[i])
        t=Time(date,format='isot',scale='tai')
        # if hasattr(t,'tai_seconds'):  #2021/11/08 HUW - this is 37 seconds too late, resort to just add constant
        #     tai=t.tai_seconds
        # else:
        tai[i]=t.unix_tai+378691163.0

    return tai

    
def anytim2cal(date_in0,tai=False,form=11,date_only=False,msec=True):

    n=np.size(date_in0)
    if isinstance(date_in0,list)==False and isinstance(date_in0,np.ndarray)==False:
        date_in=[date_in0]
    else:
        date_in=date_in0

    datemain = ["" for i in range(n)]

    for imain in range(n):

        if isinstance(date_in[imain],str) == False:
            if np.isnan(date_in[imain])==True:
                continue
            tai=True

        if tai == True: # user has supplied date_in as TAI (this is the IDL/Solarsoft TAI, or seconds since 1957/12/31 23:59:51)
            date=parse_time(date_in[imain],scale='tai',format='tai_seconds')   # convert IDL tai to date object
        else:
            date=parse_time(date_in[imain]) # assumes user has supplied date as string, in one of the accepted formats
        
        date=date.iso

        if date_only==True: #user requires only the date, without time
            # date=parse_time(date,format="iso",subformat="date") # Probably doing something dumb, but can't get output format date only like this
            date = date[0:10]
        
        if form==11:
            date = date.replace("-","/")
        elif form==8:
            date = date.replace("-","")
            date = date.replace(":","")
            date = date.replace(".","")
            date = date.replace(" ","")
            if date_only==False:
                date = date[0:14]
        else: 
            print("time_huw.anytim2cal_huw: form not recognised: ",date_in[imain],", leave unchanged")
            date = date_in[imain]

        # if user has set msec to false then check if milliseconds are in date, and remove    
        if not msec:
            isms=re.findall("[.]",date)
            if not isms==False:
                ipos = re.search("[.]", date)
                date=date[0:ipos.span()[0]]

        datemain[imain]=date

    return datemain


def yyyymmdd2cal(d0,date=False,tai=False):

    n=np.size(d0)
    if n == 1:
        d=[d0]
    else:
        d=d0

    cal = ["" for i in range(n)]
    taiout=np.empty(n)

    for i in range(n):

        isjustdate=(len(d[i]) < 9) or date==True#date2 in IDL!!!
        cal[i]=d[i][0:4]+'/'+d[i][4:6]+'/'+d[i][6:8]
        if isjustdate==False:
            cal[i]=cal[i]+' '+d[i][8:10]+':'+d[i][10:12]

        if len(d[i]) > 12:
            cal[i]=cal[i]+':'+d[i][12:14]

        if tai:
            calnow=cal[i]
            if isjustdate:
                calnow=calnow+' 00:00'
                print("time_huw.yyyymmdd2cal: you requested TAI for a date that doesn't have time, setting to 00:00")
            taiout[i]=anytim2tai(anytim2cal(calnow))

    if tai:
        calout=taiout
    else:
        calout=cal

    return cal





