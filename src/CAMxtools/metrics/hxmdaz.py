__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: hxmdaz
    :platform: Unix, Windows
    :synopsis: It reads CAMx or CMAQ output files and prepares
               grand average for the period specified.
    :details:  The sequential processes are as follows,
                1. PROCESSING COMBINE: creates indata
                2. PROCESSING MDA8O3 or MDA1O3
                3. FIND 1st or 4th highest
                4. CREATE CMAQ or CAMx output file with one time stamp:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['hxmdaz',]
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

def hxmdaz(jdbeg,jdend,outfile,comb,attr_in,nifiles,infile_h,infile_t,*,avg_hr=8,rank=4,lyyyyjjj=True,tzone=None,l2uam=False,lnew_mda8=False,tmpoutf=None):
    # Include functions
    from CAMxtools.combine.combine import combine
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools.tzone.scan_timezones import scan_timezones
    from CAMxtools.tzone.get_local_hrs import get_mda1s
    from CAMxtools.tzone.get_local_hrs import get_mda8s
    from CAMxtools.tzone.get_local_hrs import get_mda8s_a0
    from CAMxtools._cnvt._data2fin import _data2fin
    from CAMxtools.write.wrt_ioapi import wrt_ioapi
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
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

    # Set a variable name 
    spec = "O3"

    # Daily loop
    print ("  1. PROCESSING COMBINE")
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
      tracernames, indata, ovarunits = combine(None,comb,nifiles,infile,lverbose=False,loutf=False,lovarunits=True)
      s = tracernames.index(spec)
      if jdate == jdbeg : conc_hrs_utc = indata[s,:,0,:,:]
      if jdate > jdbeg  : conc_hrs_utc = np.append(conc_hrs_utc,indata[s,:,0,:,:],axis=0) #As len(s) = 1, the 2nd dimension, nt is axis=0
      jdate = int((datetime.datetime.strptime(str(jdate),"%Y%j") + datetime.timedelta(days=1)).strftime("%Y%j"))
    del indata
    print ("  2. PROCESSING MDA{} FOR ALL DAYS".format(avg_hr))
    nd = jdend - jdbeg
    if avg_hr == 1:
      mdaz = get_mda1s(conc_hrs_utc, tzone_ji, nd, nx, ny)
    elif avg_hr == 8:
      if lnew_mda8:
        mdaz = get_mda8s_a0(conc_hrs_utc, tzone_ji, nd, nx, ny)
      else:
        mdaz = get_mda8s(conc_hrs_utc, tzone_ji, nd, nx, ny)
    else:
      exit("avg_hr must be either 1 or 8.")

    # Prepare MDAZ if asked
    lout_mdaz = False
    if not tmpoutf == None: lout_mdaz = True
    if lout_mdaz:
      data2sav = np.zeros((1,nd,1,ny,nx))
      data2sav[0,:,0,:,:] = mdaz
      # Write a binary file for MDAZ
      attr_in['TSTEP']=240000
      fout = _data2fin(data2sav, tracernames, attr_in)
      if l2uam:
        wrt_uamiv(tmpoutf, fout, lsurf = True, ounits = ovarunits)
      else:
        wrt_ioapi(tmpoutf, fout, lsurf = True, ounits = ovarunits)

    # Find X highest
    if rank == 1:
      rank_name = "FIRST"
    elif rank == 4:
      rank_name = "FOURTH"
    else:
      exit("rank must be either 1 or 4.")
    print ("  3. PROCESSING {} HIGHEST".format(rank_name))
    data2sav = np.zeros((1,1,1,ny,nx))
    data2sav[0,0,0,:,:] = np.sort(mdaz,axis=0)[mdaz.shape[0]-rank,...] #data2sav[nspc,nt,nz,ny,nx]

    # Write a binary file for HXMDAZ
    fout = _data2fin(data2sav, tracernames, attr_in)
    if l2uam:
      wrt_uamiv(outfile, fout, lsurf = True, ounits = ovarunits)
    else:
      wrt_ioapi(outfile, fout, lsurf = True, ounits = ovarunits)

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
    # Check MDAZ_OUTF environment variable
    try :
      tmpoutf  = os.environ['MDAZ_OUTF']
    except :
      tmpoutf = None
    if not tmpoutf == None:
       print("  MDAZ_OUTF will be created as {}".format(tmpoutf))

    # Check arguments
    outfile  = str(sys.argv[1])
    avg_hr   = int(sys.argv[2])
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

    # Check NEW_MDA8_METHOD environment variable
    lnew_mda8 = False
    try :
      new_mda8_flag  = os.environ['NEW_MDA8_METHOD']
    except :
      new_mda8_flag  = 'F'
    if (new_mda8_flag == 'T') or (new_mda8_flag == 't') or (new_mda8_flag == 'Y') or (new_mda8_flag == 'y') :
       lnew_mda8 = True
    else:
      if avg_hr == 8:
        print("  WARNING: You ARE USING OLD MDA8 METHOD!\n 24 Running 8 hour averages will be considered to calculate MDA8 O3 instead of 17 averages")

    # Check the first day file and read
    jdate  = str(jdbeg)
    gdate  = int(datetime.datetime.strptime(str(jdate),"%Y%j").strftime("%Y%m%d"))
    if lyyyyjjj:
      infile0 = infile_h[0] + '.' + str(jdate) + '.' + infile_t[0]
    else:
      infile0 = infile_h[0] + '.' + str(gdate) + '.' + infile_t[0]
    if not os.path.exists(infile0):
      print ("{} does not exist!".format(infile0))
      print ('Program exits from hxmdaz')
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
    hxmdaz(jdbeg,jdend,outfile,comb,attr_in,nifiles,infile_h,infile_t,avg_hr=avg_hr,rank=rank,lyyyyjjj=lyyyyjjj,tzone=tz,l2uam=l2uam,lnew_mda8=lnew_mda8,tmpoutf=tmpoutf)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
