__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: run_calc_W126
    :platform: Unix, Windows
    :synopsis: It is a main script that processes W126 analysis from
               CAMx or CMAQ output files including timezone calculation.
               The processes are as follows,
                 1. Scan time zone over the domain (scan_timezones)
                 2. Calculate W126 for all timezones (calc_W126)
    :details: 1. Input files are CAMx or CMAQ binary output files
              2. Output file is IOAPI formatted file with one time stamp. 
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['run_calc_W126',]
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

def run_calc_W126(outfile,attr_in,fileh,filet,year,smo,emo,*,lyyyyjjj=True,tzone=None):

    # Include functions to call
    from CAMxtools.tzone.scan_timezones import scan_timezones
    from CAMxtools.W126.calc_W126 import calc_W126
    import numpy as np

    # Check arguments
    l1tzone = False
    if tzone != None: l1tzone = True # If users specify a specific time zone, W126 is calculated based on the time zone.

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

    # Calculate W126 per each time zone
    calc_W126(outfile,fileh,filet,year,smo,emo,tzone_ji,lyyyyjjj=lyyyyjjj,lverbose=False)


def main():
    # Include functions to call
    from CAMxtools.write.set_attr import set_attr
    import netCDF4 as ncdf4

    # Check TZONE environ variable
    try :
      tz = int(os.environ['TIMEZONE'])
    except :
      tz = None
    ## define user inputs here ##
    outfile = str(sys.argv[1])
    fileh   = str(sys.argv[2])
    filet   = str(sys.argv[3])
    year    = int(sys.argv[4])
    smo     = int(sys.argv[5])
    emo     = int(sys.argv[6])
    lyyyyjjj = False
    if str(sys.argv[7]).lower() == 'true' :
      lyyyyjjj = True # CAMx file name has yyyyjjj

    # Check the first day file and read
    jdate  = datetime.date(year,smo,1).strftime("%Y%j")
    gdate  = year*10000+smo*100+1
    if lyyyyjjj:
      infile0 = fileh + '.' + str(jdate) + '.' + filet
    else:
      infile0 = fileh + '.' + str(gdate) + '.' + filet
    if not os.path.exists(infile0):
      print ("{} does not exist!".format(infile0))
      print ('Program exits from run_calc_W126')
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

    # File attributes from an argument
    lfin0 = True
    attr_fed = {}
    attr_in = set_attr(lfin0, fin0, attr_fed)
    loutf = True

    run_calc_W126(outfile,attr_in,fileh,filet,year,smo,emo,lyyyyjjj=lyyyyjjj,tzone=tz)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
