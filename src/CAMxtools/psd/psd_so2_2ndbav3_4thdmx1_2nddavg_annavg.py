__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: psd_so2_2nddmx1_8thdmx1_2nddavg_annavg
    :platform: Unix, Windows
    :synopsis: It reads CAMx or CMAQ output files and prepares
               grand average for the period specified.
    :details:  The sequential processes are as follows,
                1. PROCESSING COMBINE: creates indata2
                2. PROCESSING METRICS FOR PREVIOUS DAY
                3. PROCESSING 2ND HIGHEST DAILY 3 HOUR BLOCK AVERAGE
                4. PROCESSING 4TH HIGHEST DAILY 1 HOUR MAXIMUM
                5. PROCESSING 2ND HIGHEST DAILY AVERAGE:
                6. PROCESSING GRAND AVERAGE:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['psd_so2_2nddavg_annavg',]
import sys
if sys.version_info.major == 3:
    from io import BytesIO as StringIO
    commaspace = u', '
    semicolon = b';'
else:
    from StringIO import StringIO
    commaspace = ', '
    semicolon = ';'
    BrokenPipeError = None
from pathlib import Path
import os
import datetime

def psd_so2_2ndbav3_4thdmx1_2nddavg_annavg(jdbeg,jdend,out_anndavg,out_2nddavg,out_2ndbav3,out_4thdmx1,csv_anndavg,csv_2nddavg,csv_2ndbav3,csv_4thdmx1,comb,attr_in,nifiles,infile_h,infile_t,*,lyyyyjjj=True,tzone=None,l2uam=False,tmpoutf=None):
    # Include functions
    from CAMxtools.combine.combine import combine
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools.tzone.scan_timezones import scan_timezones
    from CAMxtools.tzone.get_local_hrs import get_davg
    from CAMxtools.tzone.get_local_hrs import get_dbav3
    from CAMxtools.tzone.get_local_hrs import get_mda1
    from CAMxtools._cnvt._data2fin import _data2fin
    from CAMxtools.write.wrt_ioapi import wrt_ioapi
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.psd.wrt_csv_for_psd import wrt_csv_for_psd
    import numpy as np

    # Check arguments
    l1tzone = False
    if tzone != None: l1tzone = True # If users specify a specific time zone, MATS MDA8 O3 is calculated based on the time zone.

    # Scan time zone over the domain
    tzfile = None
    dum, tzone_stlji = scan_timezones(tzfile, attr_in, loutf = False)
    tzone_ji = tzone_stlji[0,0,0,:,:].astype(np.int)
    tzones = np.unique(tzone_ji)
    if l1tzone:
      print("A SINGLE TIMEZONE, {} will be applied".format(tzone))
      if tzone in tzones:
        ny = attr_in['NROWS']
        nx = attr_in['NCOLS']
        tzone_ji = np.zeros((ny,nx)).astype(np.int) + tzone
        tzones = [tzone]
      else:
        exit("YOUR TIMEZONE SPECIFIED IS OUT OF DOMAIN")

    # Get attribute
    nx  = attr_in['NCOLS']
    ny  = attr_in['NROWS']

    # Daily loop
    jdate = jdbeg
    while (jdate <= jdend):
      print ("Processing {}".format(jdate))
      gdate = int(datetime.datetime.strptime(str(jdate),"%Y%j").strftime("%Y%m%d"))
      infile = []
      for ifile in range(0,nifiles):
        if lyyyyjjj:
          infile.append(infile_h[ifile]+'.'+str(jdate)+'.'+infile_t[ifile])
        else:
          infile.append(infile_h[ifile]+'.'+str(gdate)+'.'+infile_t[ifile])
      print ("  1. PROCESSING COMBINE")
      tracernames, indata2, ovarunits = combine(None,comb,nifiles,infile,lverbose=False,loutf=False,lovarunits=True)
      if (jdate > jdbeg):
        print ("  2. PROCESSING METRICS FOR PREVIOUS DAY")
        concs_48_utc = np.append(indata1[:,:,0,:,:],indata2[:,:,0,:,:],axis=1)
        davg_1day = get_davg(concs_48_utc, tzone_ji, nx, ny)
        dbav3_1day = get_dbav3(concs_48_utc, tzone_ji, nx, ny) # dbav3_1day[8,nspc,ny,nx]
        dmda1_1day, dind_hr_1day = get_mda1(concs_48_utc, tzone_ji, nx, ny) # dmda1_1day[nspc,ny,nx], dind_hr_1day[nspc,ny,nx]
        if jdate == jdbeg + 1:
          davgs = np.array([davg_1day]) #davgs[nd,nspc,ny,nx]
          dbav3s = dbav3_1day #dbav3s[nd*8,nspc,ny,nx]
          dmda1s = np.array([dmda1_1day]) #dmda1s[nd,nspc,ny,nx]
          dind_hrs = np.array([dind_hr_1day]) #dind_hrs[nd,nspc,ny,nx]
        if jdate > jdbeg + 1:
          davgs = np.append(davgs,np.array([davg_1day]),axis=0)
          dbav3s = np.append(dbav3s,dbav3_1day,axis=0)
          dmda1s = np.append(dmda1s,np.array([dmda1_1day]),axis=0)
          dind_hrs = np.append(dind_hrs,np.array([dind_hr_1day]),axis=0)
      indata1 = indata2
      jdate = int((datetime.datetime.strptime(str(jdate),"%Y%j") + datetime.timedelta(days=1)).strftime("%Y%j"))
    del indata2
    del indata1

    # Set nspc
    nspc = len(tracernames)

    # Calculate 2nd highest daily 8 hour block average
    print ("  3. PROCESSING 2ND HIGHEST 3 HOUR BLOCK AVERAGE")
    data2sav = np.zeros((nspc,1,1,ny,nx))
    rank = 2
    indX8 = np.argsort(dbav3s,axis=0)[dbav3s.shape[0]-rank,...] #data2sav[nspc,nt,nz,ny,nx], dbav3s[nd*8,nspc,ny,nx], indX8[nsp,ny,nx]
    gr = np.ogrid[0:dbav3s.shape[0],0:dbav3s.shape[1],0:dbav3s.shape[2],0:dbav3s.shape[3]]
    gr[0] = indX8
    data2sav[:,0,0,:,:] = dbav3s[gr]
    # Write a csv file
    jdays = attr_in['SDATE']+indX8[0,:,:]//8
    hours = np.zeros((ny,nx))+(indX8[0,:,:]%8)*3 # 3 hour block average
    data2csv = np.zeros((nspc,ny,nx))
    data2csv = data2sav[:,0,0,:,:]
    wrt_csv_for_psd( csv_2ndbav3, tracernames, data2csv, jdays, hours)
    # Write a binary file
    fout = _data2fin(data2sav, tracernames, attr_in)
    if l2uam:
      wrt_uamiv(out_2ndbav3, fout, lsurf = True, ounits = ovarunits)
    else:
      wrt_ioapi(out_2ndbav3, fout, lsurf = True, ounits = ovarunits)

    # Calculate 4th highest daily 1 hour maximum
    print ("  4. PROCESSING 4TH HIGHEST DAILY 1 HOUR MAXIMUM")
    data2sav = np.zeros((nspc,1,1,ny,nx))
    rank = 4
    ind = np.argsort(dmda1s,axis=0)[dmda1s.shape[0]-rank,...] #data2sav[nspc,nt,nz,ny,nx], dmda1s[nd,nspc,ny,nx], ind[nsp,ny,nx]
    gr = np.ogrid[0:dmda1s.shape[0],0:dmda1s.shape[1],0:dmda1s.shape[2],0:dmda1s.shape[3]]
    gr[0] = ind
    data2sav[:,0,0,:,:] = dmda1s[gr] #dmda1s[gr].shape = (1,nspc,ny,nx)
    # Write a csv file
    jdays = attr_in['SDATE']+ind[0,:,:]
    hours = np.zeros((ny,nx))+dind_hrs[gr][0,0,:,:] #dind_hrs[gr].shape = (1,nspc,ny,nx)
    data2csv = np.zeros((nspc,ny,nx))
    data2csv = data2sav[:,0,0,:,:]
    wrt_csv_for_psd( csv_4thdmx1, tracernames, data2csv, jdays, hours)
    # Write a binary file
    fout = _data2fin(data2sav, tracernames, attr_in)
    if l2uam:
      wrt_uamiv(out_4thdmx1, fout, lsurf = True, ounits = ovarunits)
    else:
      wrt_ioapi(out_4thdmx1, fout, lsurf = True, ounits = ovarunits)

    # Calculate 2nd highest daily average
    print ("  5. PROCESSING 2ND HIGHEST DAILY AVERAGE")
    data2sav = np.zeros((nspc,1,1,ny,nx))
    rank = 2
    ind = np.argsort(davgs,axis=0)[davgs.shape[0]-rank,...] #data2sav[nspc,nt,nz,ny,nx], davgs[nd,nspc,ny,nx], ind[nsp,ny,nx]
    gr = np.ogrid[0:davgs.shape[0],0:davgs.shape[1],0:davgs.shape[2],0:davgs.shape[3]]
    gr[0] = ind
    data2sav[:,0,0,:,:] = davgs[gr]
    # Write a csv file
    jdays = attr_in['SDATE']+ind[0,:,:] # Annual average, set the beginning date 
    hours = np.zeros((ny,nx)) # daily average, set zeros
    data2csv = np.zeros((nspc,ny,nx))
    data2csv = data2sav[:,0,0,:,:]
    wrt_csv_for_psd( csv_2nddavg, tracernames, data2csv, jdays, hours)
    # Write a binary file
    fout = _data2fin(data2sav, tracernames, attr_in)
    if l2uam:
      wrt_uamiv(out_2nddavg, fout, lsurf = True, ounits = ovarunits)
    else:
      wrt_ioapi(out_2nddavg, fout, lsurf = True, ounits = ovarunits)

    # Calculate grand average
    print ("  6. PROCESSING GRAND AVERAGE")
    data2sav = np.zeros((nspc,1,1,ny,nx))
    data2sav[:,0,0,:,:] = np.mean(davgs,axis=0) #data2sav[nspc,nt,nz,ny,nx], davgs[nd,nspc,ny,nx]
    # Write a csv file
    jdays = np.zeros((ny,nx))+attr_in['SDATE'] # Annual average, set the beginning date 
    hours = np.zeros((ny,nx)) # Annual average, set zeros
    data2csv = np.zeros((nspc,ny,nx))
    data2csv = data2sav[:,0,0,:,:]
    wrt_csv_for_psd( csv_anndavg, tracernames, data2csv, jdays, hours)
    # Write a binary file
    fout = _data2fin(data2sav, tracernames, attr_in)
    if l2uam:
      wrt_uamiv(out_anndavg, fout, lsurf = True, ounits = ovarunits)
    else:
      wrt_ioapi(out_anndavg, fout, lsurf = True, ounits = ovarunits)

def main():
    # Include functions to call
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr

    # Check TZONE environment variable
    try :
      tz = int(os.environ['TIMEZONE'])
    except :
      tz = None
    # Check OUT2UAM environment variable
    l2uam = False
    try :
      uamflag  = os.environ['OUT2UAM']
    except :
      uamflag  = 'F'
    if (uamflag == 'T') or (uamflag == 't') or (uamflag == 'Y') or (uamflag == 'y') :
       l2uam = True

    # Check arguments
    out_anndavg = str(sys.argv[1])
    out_2nddavg = str(sys.argv[2])
    out_2ndbav3 = str(sys.argv[3])
    out_4thdmx1 = str(sys.argv[4])
    csv_anndavg = str(sys.argv[5])
    csv_2nddavg = str(sys.argv[6])
    csv_2ndbav3 = str(sys.argv[7])
    csv_4thdmx1 = str(sys.argv[8])
    jdbeg = int(sys.argv[9])
    jdend = int(sys.argv[10])
    lyyyyjjj = False
    if str(sys.argv[11]).lower() == "true" :
      lyyyyjjj = True # CAMx file name has yyyyjjj
    comb = str(sys.argv[12])
    nifiles = int(sys.argv[13])
    infile_h = []
    infile_t = []
    # Build input file header and tail lists
    nargs2drop = 14 # No. of arguments up to "infile1 header" (inclusive)
    for ifile in range(0,nifiles):
      infile_h.append(str(sys.argv[nargs2drop+(2*ifile)]))
      infile_t.append(str(sys.argv[nargs2drop+1+(2*ifile)]))

    # Check the first day file and read
    jdate  = str(jdbeg)
    gdate  = int(datetime.datetime.strptime(str(jdate),"%Y%j").strftime("%Y%m%d"))
    if lyyyyjjj:
      infile0 = infile_h[0] + '.' + str(jdate) + '.' + infile_t[0]
    else:
      infile0 = infile_h[0] + '.' + str(gdate) + '.' + infile_t[0]
    if not os.path.exists(infile0):
      print ("{} does not exist!".format(infile0))
      print ('Program exits from psd_so2_2nddavg_annavg')
      exit()
    try:
      fin0 = uamiv(infile0)
    except:
      try:
        fin0 = ncdf4.Dataset(infile0)
      except:
        print ("Unrecognized file type")
        print ("infile1 = {}".format(infile0))
        exit()
    
    # Set file attributes
    lfin0 = True
    attr_fed = {}
    attr_in = set_attr(lfin0,fin0,attr_fed)

    # Do the main process
    psd_so2_2ndbav3_4thdmx1_2nddavg_annavg(jdbeg,jdend,out_anndavg,out_2nddavg,out_2ndbav3,out_4thdmx1,csv_anndavg,csv_2nddavg,csv_2ndbav3,csv_4thdmx1,comb,attr_in,nifiles,infile_h,infile_t,lyyyyjjj=lyyyyjjj,tzone=tz,l2uam=l2uam)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
