from sunpy.time import parse_time
from astropy.time import Time
import numpy as np

def anytim2tai_huw(date_in):

    date=parse_time(date_in)
    t=Time(date,format='isot',scale='tai')
    tai=t.tai_seconds

    return tai

    
def anytim2cal_huw(date_in,tai=False,form=11,date_only=False):

    if isinstance(date_in,str) == False:
        tai=True

    if tai == True: # user has supplied date_in as TAI (this is the IDL/Solarsoft TAI, or seconds since 1957/12/31 23:59:51)
        date=parse_time(date_in,scale='tai',format='tai_seconds')   # convert IDL tai to date object
    else:
        date=parse_time(date_in) # assumes user has supplied date as string, in one of the accepted formats
    
    date=date.iso

    if date_only==True: #user requires only the date, without time
        # date=parse_time(date,format="iso",subformat="date") # Probably doing something dumb, but can't get output format date only like this
        date = [date[i][0:10] for i in range(np.size(date))]
    
    if form==11:
        date = [date[i].replace("-","/") for i in range(np.size(date))]
    elif form==8:
        date = [date[i].replace("-","") for i in range(np.size(date))]
        date = [date[i].replace(":","") for i in range(np.size(date))]
        date = [date[i].replace(".","") for i in range(np.size(date))]
        date = [date[i].replace(" ","") for i in range(np.size(date))]
        if date_only==False:
            date = [date[i][0:14] for i in range(np.size(date))]
    else: 
        print("time_huw.anytim2cal_huw: form not recognised")
        return

    return date
    




