
def phys_constants(system='cgs',rsun=False,au=False):

    if system == 'cgs':
        c={ 
            "Rsun":6.955e+10    ,
            "AU": 1.495978707e+13,
            "boltzK":1.380622e-16,
            "clight":2.997925e+10,
            "hPlanck":6.62607004e-27,
            "eV":1.6022e-12,
            "mprot":1.6726219e-24,
            "melect":9.10938356e-28,
            "eelect":4.8032e-10, 
            "gGrav":6.67408e-8,
            "Msun":1.98892e33
        }
    elif system == 'si':
        c={	
            "Rsun":6.955e+8		, 
            "AU": 1.495978707e+11, 
            "boltzK":1.380622e-23	, 
            "clight":2.997925e+8	, 
            "hPlanck":6.62607004e-34	, 
            "eV":1.6022e-19		, 
            "mprot":1.6726219e-27	, 
            "melect":9.10938356e-31	, 
            "eelect":1.60217662e-19 	, 
            "gGrav":6.67408e-11		, 
            "Msun":1.98892e30	
        } 				


    if rsun==True:
        c=c["Rsun"]
    elif au==True:
        c=c["AU"]

    return c
