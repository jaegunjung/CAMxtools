__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: prep_pm_mats
    :platform: Unix, Windows
    :synopsis: It reads CAMx or CMAQ output files and prepares
               O3 MATS input file, which is a csv file with
               MDA8 values at each grid cells.
    :details:  The sequential processes are as follows,
                1. PROCESSING COMBINE: creates indata2
                2. PROCESSING MDA8:
                3. CREATE output file:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['prep_pm_mats',]
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

def prep_pm_mats(jdbeg,jdend,outfile,comb,attr_in,nifiles,infile_h,infile_t,*,lyyyyjjj=True,tzone=None,lno_edges=True,convfac=1.):
    # Include functions
    from CAMxtools.combine.combine import combine
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools.tzone.scan_timezones import get_lcc
    from CAMxtools.regrid.projection import ll2latlon
    from CAMxtools.regrid.projection import lcp2latlon
    from CAMxtools.tzone.scan_timezones import scan_timezones
    from CAMxtools.tzone.get_local_hrs import get_davg
    from CAMxtools.MATS.wrt_csv_for_pm_mats import wrt_csv_for_pm_mats
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
    dxy = attr_in['XCELL']
    nx  = attr_in['NCOLS']
    ny  = attr_in['NROWS']
    x0 = attr_in['XORIG']
    y0 = attr_in['YORIG']

    # Get lat and lon at the center of grid cells
    if attr_in['GDTYP'] == 1:
      lats, lons = ll2latlon (x0, y0, dxy, nx, ny, lno_edges = lno_edges)
    elif attr_in['GDTYP'] == 2:
      lcc = get_lcc(attr_in)
      lats, lons = lcp2latlon (lcc, dxy, nx, ny, lno_edges = lno_edges)
    else:
      print("Your GDTYP is {}".format(attr_in['GDTYP']))
      exit("This program currently supports LATLON (GDTYP = 1) and LCP (GDTYP = 2) only.")

    # Daily loop
    jdays = []
    jdate = jdbeg
    while (jdate <= jdend):
      jdays.append(jdate)
      print ("Processing {}".format(jdate))
      gdate = int(datetime.datetime.strptime(str(jdate),"%Y%j").strftime("%Y%m%d"))
      infile = []
      for ifile in range(0,nifiles):
        if lyyyyjjj:
          infile.append(infile_h[ifile]+'.'+str(jdate)+'.'+infile_t[ifile])
        else:
          infile.append(infile_h[ifile]+'.'+str(gdate)+'.'+infile_t[ifile])
      print ("  1. PROCESSING COMBINE")
      tracernames, indata2 = combine(None,comb,nifiles,infile,lverbose=False,loutf=False)
      if (jdate > jdbeg):
        print ("  2. PROCESSING DAILY AVERAGE FOR PREVIOUS DAY")
        concs_48_utc = np.append(indata1[:,:,0,:,:],indata2[:,:,0,:,:],axis=1)
        davg_1day = get_davg(concs_48_utc, tzone_ji, nx, ny)
        if jdate == jdbeg + 1: davgs = np.array([davg_1day]) #davgs[nd,nspc,ny,nx]
        if jdate > jdbeg + 1: davgs = np.append(davgs,np.array([davg_1day]),axis=0)
      indata1 = indata2
      jdate = int((datetime.datetime.strptime(str(jdate),"%Y%j") + datetime.timedelta(days=1)).strftime("%Y%j"))
    del indata2
    del indata1
    
    # Exclude buffer cells if needed
    if lno_edges:
      davgs_chked  = davgs[:,:,1:-1,1:-1]
    else:
      davgs_chked  = davgs

    # Write csv file
    wrt_csv_for_pm_mats(outfile,davgs_chked,jdays,lats,lons,tracernames,convfac)

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
    # Check INCLUDE_BUFFER_CELLS environment variable
    lno_edges = True
    try:
      bfflag = os.environ['INCLUDE_BUFFER_CELLS']
    except:
      bfflag = 'F'
    if (bfflag.lower() == 't') or (bfflag.lower() == 'y'):
      lno_edges = False
    # Check CONVFAC environment variable
    try :
      convfac = float(os.environ['CONVFAC'])
    except :
      convfac = 1.

    # Check arguments
    outfile = str(sys.argv[1])
    jdbeg = int(sys.argv[2])
    jdend = int(sys.argv[3])
    lyyyyjjj = False
    if str(sys.argv[4]).lower() == "true" :
      lyyyyjjj = True # CAMx file name has yyyyjjj
    comb = str(sys.argv[5])
    nifiles = int(sys.argv[6])
    infile_h = []
    infile_t = []
    # Build input file header and tail lists
    nargs2drop = 7 # No. of arguments up to "infile1 header" (inclusive)
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
      print ('Program exits from prep_pm_mats')
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

    # Check whether buffer cell is 0. when lno_edge is off.
    spc0 = getattr(fin0,'VAR-LIST').split()[0]
    val0 = fin0.variables[spc0][0,0,0,0]
    if lno_edges:
      if val0 != 0.0:
        print("WARNING: The value at the SW corner cell is not zero while you want to exclude buffer cells.")
    else:
      assert val0 != 0.0, "The value at the SW corner cell is zero while you want to include buffer cells. Program exits."

    # Set file attributes
    lfin0 = True
    attr_fed = {}
    attr_in = set_attr(lfin0,fin0,attr_fed)

    # Do the main process
    prep_pm_mats(jdbeg,jdend,outfile,comb,attr_in,nifiles,infile_h,infile_t,lyyyyjjj=lyyyyjjj,tzone=tz,lno_edges=lno_edges,convfac=convfac)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
