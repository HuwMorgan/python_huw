from sunpy.time import parse_time
from astropy.time import Time
import numpy as np

def anytim2tai(date_in):

    date=parse_time(date_in)
    t=Time(date,format='isot',scale='tai')
    if hasattr(t,'tai_seconds'):
        tai=t.tai_seconds
    else:
        tai=t.unix_tai+378691200.0

    return tai

    
def anytim2cal(date_in0,tai=False,form=11,date_only=False):

    n=np.size(date_in0)
    if n == 1:
        date_in=[date_in0]
    else:
        date_in=date_in0

    datemain=[]

    for imain in range(n):

        if isinstance(date_in[imain],str) == False:
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

        datemain.append(date)

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





