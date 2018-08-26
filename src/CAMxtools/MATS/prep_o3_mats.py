__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: prep_o3_mats
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

__all__=['prep_o3_mats',]
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

def prep_o3_mats(jdbeg,jdend,outfile,comb,attr_in,nifiles,infile_h,infile_t,*,lyyyyjjj=True,tzone=None,lno_edges=True,convfac=1000.,lnew_mda8=False):
    # Include functions
    from CAMxtools.combine.combine import combine
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools.tzone.scan_timezones import get_lcc
    from CAMxtools.regrid.projection import ll2latlon
    from CAMxtools.regrid.projection import lcp2latlon
    from CAMxtools.tzone.scan_timezones import scan_timezones
    from CAMxtools.tzone.get_local_hrs import get_mda8s
    from CAMxtools.tzone.get_local_hrs import get_mda8s_a0
    from CAMxtools.MATS.wrt_csv_for_o3_mats import wrt_csv_for_o3_mats
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

    # Set a variable name
    spec = "O3"

    # Daily loop
    print ("  1. PROCESSING COMBINE")
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
      tracernames, indata = combine(None,comb,nifiles,infile,lverbose=False,loutf=False)
      s = tracernames.index(spec)
      if jdate == jdbeg : conc_hrs_utc = indata[s,:,0,:,:]
      if jdate > jdbeg  : conc_hrs_utc = np.append(conc_hrs_utc,indata[s,:,0,:,:],axis=0) #As len(s) = 1, the 2nd dimension, nt is axis=0
      jdate = int((datetime.datetime.strptime(str(jdate),"%Y%j") + datetime.timedelta(days=1)).strftime("%Y%j"))
    del indata
    print ("  2. PROCESSING MDA8 FOR ALL DAYS")
    nd = jdend - jdbeg
    if lnew_mda8:
      mda8 = get_mda8s_a0(conc_hrs_utc, tzone_ji, nd, nx, ny)
    else:
      mda8 = get_mda8s(conc_hrs_utc, tzone_ji, nd, nx, ny)
    
    # Exclude buffer cells if needed
    if lno_edges:
      mda8_chked  = mda8[:,1:-1,1:-1]
    else:
      mda8_chked  = mda8

    # Write csv file
    wrt_csv_for_o3_mats(outfile,mda8_chked,jdays,lats,lons,convfac)

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
    # Check NEW_MDA8_METHOD environment variable
    lnew_mda8 = False
    try :
      new_mda8_flag  = os.environ['NEW_MDA8_METHOD']
    except :
      new_mda8_flag  = 'F'
    if (new_mda8_flag == 'T') or (new_mda8_flag == 't') or (new_mda8_flag == 'Y') or (new_mda8_flag == 'y') :
       lnew_mda8 = True
    else:
       print("  WARNING: You ARE USING OLD MDA8 METHOD!\n 24 Running 8 hour averages will be considered to calculate MDA8 O3 instead of 17 averages")
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
      convfac = 1000.

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
      print ('Program exits from prep_o3_mats')
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
    prep_o3_mats(jdbeg,jdend,outfile,comb,attr_in,nifiles,infile_h,infile_t,lyyyyjjj=lyyyyjjj,tzone=tz,lno_edges=lno_edges,convfac=convfac,lnew_mda8=lnew_mda8)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
