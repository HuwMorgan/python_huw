def timerange2datacube(middate,mintrange=None, maxtimegap=None):
    
    """ Process the first step of tomography
        Parameters
        ----------
        startdate : `str`
        enddate : `str`
        mintrange : number of seconds
        maxtimegap : number of seconds
        middate : `str`
        synthetic : `bool`, optional default `False`
        usrdir : `str`, optional
        save_data : `bool`, optional default `True`
        Returns
        -------
            a datacube
        Details
        -------
        - open ~two weeks of pB fits files and save in a datacube
        - Extract one ’slice’ of the data at a constant distance from the Sun
        - Prepare some other variables and arrays that are needed
          for the tomography (e.g. set up lines of sight and
          Carrington co-ordinates etc).
        Written by: Huw Morgan
    """

    SEC_IN_HOUR = 3600
    MAX_N_FILES = 50

    startdate,enddate = return_timerange(middate)
    
    print(startdate, '=>', enddate)

    
    stereo='a'
    instr = 'cor2'
    spacecraft = stereo
    mission = 'stereo'
    # topdir = os.getenv('PROCESSED_DATA')+'/stereo/secchi'
    # caldir = topdir+'/separated_new/'+spacecraft+'/cor2/dat'
    # srchstr = instr+'_'+spacecraft+'_stereo_quiescnt_*.dat'
    topdir = os.path.join(os.getenv('PROCESSED_DATA'),
                            'stereo/secchi/pb', spacecraft)
    caldir = os.path.join(os.getenv('PROCESSED_DATA'),
                            'stereo/secchi/pb', spacecraft, 'cor2')
    # srchstr = instr + '_' + spacecraft + '_bk_*.dat'
    srchstr = 'cor2_a_pb_*.pkl'
    system = 'sta'

    starttai = anytim2tai(startdate)
    endtai = anytim2tai(enddate)
    range_ = endtai - starttai

    dates=timegrid(startdate,enddate,delta=1,days=True)
    days=anytim2cal(dates,form=11,date_only=True)
    
    ndays = len(days)


    cntfiles = 0
    files = np.empty(0,dtype=object)
    for iday in range(ndays):
        fnow = glob.glob(os.path.join(caldir, days[iday],'')+srchstr)
        cnt = len(fnow)
        if cnt == 0:
            continue
        else:
            files=np.append(files,fnow)
            cntfiles += cnt

    if cntfiles == 0:
        logging.error('No files found (timerange2datacube)')
        return -1

    if cntfiles < MAX_N_FILES:
        logging.error('Not enough files found (timerange2datacube)')
        return -1

    np.sort(files)

    for ifile in range(cntfiles):
        if ifile % 50 == 0:
            print(files[ifile])
        
        afile=open(files[ifile],"rb")
        d=pickle.load(afile)
        afile.close()

        # d = json.load(files[ifile])

        if ifile == 0:
            imm = np.zeros((d['npa'], d['nht'], cntfiles))
            dates = np.empty(cntfiles,dtype=object)
            filesorig = np.empty(cntfiles,dtype=object)

        imm[:, :, ifile] = np.transpose(d['pb'])
        dates[ifile]=d['date']
        filesorig[ifile]=d['files'][0]

    t = np.sum(np.sum(imm, axis=0), axis=0)
    m = ndimage.median_filter(t, size=15)
    df = np.abs(t - m)/m
    cntfilesorig = cntfiles
    mad = median_absolute_deviation(df, nan_policy='omit')
    index=df < mad*7
    indok = np.squeeze(np.where(index))
    cntfiles = np.count_nonzero(index)


    if cntfiles == 0:
        logging.error('Very strange problem with data!')
        return

    print('Number of rejected files =', cntfilesorig - cntfiles)
    imm = imm[:, :, indok]
    dates = dates[indok]
    files = files[indok]
    filesorig = filesorig[indok]

    if mintrange is not None:
        dtai = anytim2tai(dates)

        if (max(dtai) - min(dtai)) < mintrange:
            logging.error('Available data time range'
                        'is too small (timerange2datacube)')
            return

        if maxtimegap is not None:
            dt = dtai[1:] - dtai[0:-1]
            cntbad = np.count_nonzero(dt >= maxtimegap)
            if cntbad > 0:
                logging.error('There are data gaps larger than maxtimegap =',
                            maxtimegap/SEC_IN_HOUR, 'hours (timerange2datacube)')
                return 

    imm[0,3:5,10]=-1
    indexbad=np.logical_or(imm < 0,~np.isfinite(imm))
    cntbad = np.count_nonzero(indexbad)
    if cntbad > 0:
        imm = np.where(indexbad,np.nan,imm)
        print(cntbad, 'bad pixels (all set now to NAN)')


    d = {'dates': list(dates), 'files': list(files), 'im': imm, 'nr': d['nht'],
         'npa': d['npa'], 'n': cntfiles, 'rra': d['rra'], 'para': d['para'],
         'filesorig': list(filesorig), 'geometry': 'polar', 'system': "sta",
         'type': 'pb'}

    return d