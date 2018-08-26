__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: hxdavg_naaqs
    :platform: Unix, Windows
    :synopsis: It reads CAMx or CMAQ output files and prepares
               grand average for the period specified.
    :details:  The sequential processes are as follows,
                1. PROCESSING COMBINE: creates indata2
                2. PROCESSING DAILY AVERAGE:
                3. FIND 1st or 8th highest or grand average
                4. CREATE CMAQ or CAMx output file with one time stamp:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['hxdavg_naaqs',]
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

def hxdavg_naaqs(jdbeg,jdend,outfile,csvfile,comb,attr_in,nifiles,infile_h,infile_t,*,rank=8,lyyyyjjj=True,tzone=None,l2uam=False,tmpoutf=None):
    # Include functions
    from CAMxtools.combine.combine import combine
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools.tzone.scan_timezones import scan_timezones
    from CAMxtools.tzone.get_local_hrs import get_davg
    from CAMxtools._cnvt._data2fin import _data2fin
    from CAMxtools.write.wrt_ioapi import wrt_ioapi
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.naaqs.wrt_csv_for_naaqs import wrt_csv_for_naaqs
    import numpy as np
    import gc

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

    # Set a variable name 
    spec = "PM25"

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
        # One time execution - Check 'O3' is the first output species and index output trcnames which excludes O3 from tracernames
        if jdate == jdbeg + 1:
          s = tracernames.index(spec)
          if s != 0: exit('PM25 must be the first species in combine spec_def')
          ntracers = len(tracernames)
          trcnames = tracernames[1:]
        print ("  2. PROCESSING METRICS FOR PREVIOUS DAY")
        concs_48_utc = np.append(indata1[:,:,0,:,:],indata2[:,:,0,:,:],axis=1)
        davg_1day = get_davg(concs_48_utc, tzone_ji, nx, ny)
        if jdate == jdbeg + 1: np.array([davg_1day]) #davgs[nd,nspc,ny,nx]
        if jdate > jdbeg + 1: davgs = np.append(davgs,np.array([davg_1day]),axis=0)
      indata1 = indata2
      jdate = int((datetime.datetime.strptime(str(jdate),"%Y%j") + datetime.timedelta(days=1)).strftime("%Y%j"))
      gc.collect()
    del indata2
    del indata1
    gc.collect()

    # Prepare MAVG if asked
    nspc = len(tracernames)
    lout_davg = False
    if not tmpoutf == None: lout_davg = True
    if lout_davg:
      nd = jdend - jdbeg
      data2sav = np.zeros((ntracers,nd,1,ny,nx))
      data2sav[:,:,0,:,:] = np.einsum('jikl->ijkl',davgs)
      # Write a binary file for DAVG
      attr_in['TSTEP']=240000
      fout = _data2fin(data2sav, tracernames, attr_in)
      if l2uam:
        wrt_uamiv(tmpoutf, fout, lsurf = True, ounits = ovarunits)
      else:
        wrt_ioapi(tmpoutf, fout, lsurf = True, ounits = ovarunits)
      gc.collect()

    # Find X highest
    if rank == 1:
      rank_name = "FIRST"
    elif rank == 8:
      rank_name = "EIGHTH"
    elif rank == 0:
      rank_name = "AVERAGE"
    else:
      exit("rank must be either 1 or 8 or 0.")
    if rank == 0:
      print ("  3. PROCESSING AVERAGE")
    else:
      print ("  3. PROCESSING {} HIGHEST".format(rank_name))
    data2sav = np.zeros((ntracers,1,1,ny,nx))
    nd = jdend - jdbeg
    if (nd < rank):
      print('Your no. of days is less than the rank you specified.')
      print('No. of days from input files = {}'.format(nd))
      print('Rank you select is = {}'.format(rank))
      exit('Either reduce your rank or increase no. of days from input files')
    if rank == 0:
      data2sav[:,0,0,:,:] = np.mean(davgs,axis=0) #data2sav[nspc,nt,nz,ny,nx], davgs[nd,nspc,ny,nx]
    else:
      davg_pm25 = davgs[:,0,:,:] #davg_pm25[nd,ny,nx], davgs[nd,nspc,ny,nx]
      ind = np.argsort(davg_pm25,axis=0)[nd-rank,...] #davgs[nd,nspc,ny,nx], ind[ny,nx]
      gr = np.ogrid[0:davg_pm25.shape[0],0:davg_pm25.shape[1],0:davg_pm25.shape[2]]
      gr[0] = ind
      data2sav[0,0,0,:,:] = davg_pm25[gr] #data2sav[nspc,nt,nz,ny,nx]
      for i in range(nx):
        for j in range(ny):
          data2sav[1:,0,0,j,i] = davgs[ind[j,i],1:,j,i]
    # Write a csv file
    hours = np.zeros((ny,nx)) # Average of DAVG, set zeros.
    if rank == 0:
      jdays = np.zeros((ny,nx))+attr_in['SDATE'] # Average of DAVG, set the beginning date
    else:
      jdays = attr_in['SDATE']+ind[:,:]
    data2csv = np.zeros((ntracers,ny,nx))
    data2csv = data2sav[:,0,0,:,:]
    wrt_csv_for_naaqs( csvfile, tracernames, data2csv, jdays, hours)
    gc.collect()
    # Write a binary file for HXDAVG
    fout = _data2fin(data2sav, tracernames, attr_in)
    if l2uam:
      wrt_uamiv(outfile, fout, lsurf = True, ounits = ovarunits)
    else:
      wrt_ioapi(outfile, fout, lsurf = True, ounits = ovarunits)
    gc.collect()

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
    # Check DAVG_OUTF environment variable
    try :
      tmpoutf  = os.environ['DAVG_OUTF']
    except :
      tmpoutf = None
    if not tmpoutf == None:
       print("  DAVG_OUTF will be created as {}".format(tmpoutf))

    # Check arguments
    outfile  = str(sys.argv[1])
    csvfile  = str(sys.argv[2])
    rank     = int(sys.argv[3])
    jdbeg    = int(sys.argv[4])
    jdend    = int(sys.argv[5])
    lyyyyjjj = False
    if str(sys.argv[6]).lower() == "true" :
      lyyyyjjj = True # CAMx file name has yyyyjjj
    comb     = str(sys.argv[7])
    nifiles  = int(sys.argv[8])
    infile_h = []
    infile_t = []
    # Build input file header and tail lists
    nargs2drop = 9 # No. of arguments up to "infile1 header" (inclusive)
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
      print ('Program exits from hxdavg_naaqs')
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
    hxdavg_naaqs(jdbeg,jdend,outfile,csvfile,comb,attr_in,nifiles,infile_h,infile_t,rank=rank,lyyyyjjj=lyyyyjjj,tzone=tz,l2uam=l2uam,tmpoutf=tmpoutf)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
