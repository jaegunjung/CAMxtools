_ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: scan_timezones
    :platform: Unix, Windows
    :synopsis: Takes a sample IOAPI file for the domain definition and
               output file attributes. It create an output IOAPI file
               which shows timezone in each grid cell.
    :warning:  tzwhere and pytz are not in Anaconda3. Install it using
               a guidance desribed in ../doc/0dependency.txt.
    :history:

... moduleauthor:: Jaegun Jung <jjung@ramboll.com>
"""

__all__=['scan_timezones','proj_latlon_single','proj_ij_single','get_lcc']
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
from PseudoNetCDF.camxfiles.Memmaps import uamiv
from pyproj import Proj
import calendar
import numpy as np
import netCDF4 as ncdf4
import math
import datetime
import os

def proj_latlon_single (i, j, dxy, lcc):
    xcoord = np.float32(i) * dxy + dxy/2 ; ycoord = np.float32(j) * dxy + dxy/2
    lon,lat = lcc(xcoord,ycoord,inverse = 'true')
    return lat, lon

def proj_ij_single (lon, lat, dxy, lcc):
    lcpx0,lcpy0 = lcc(lon,lat) # lcpx = lcpx0 + x0, lcpy = lcpy0 + y0
    i = int(lcpx0/dxy) + 1
    j = int(lcpy0/dxy) + 1
    return i, j

def get_lcc (attr_in):
    # Extract lambert projection parameters and dxy to convert to lat/lon
    lon0 = str(attr_in['XCENT'])
    lat0 = str(attr_in['YCENT'])
    lat1 = str(attr_in['P_ALP'])
    lat2 = str(attr_in['P_BET'])
    # x0, y0 need to have opposite signs to the original definition to work with the Proj package.
    x0   = str(-attr_in['XORIG'])
    y0   = str(-attr_in['YORIG'])
    #global lcc 
    lcc = Proj('+proj=lcc +a=6370000, +b=6370000, +lon_0='+lon0+' +lat_0='+lat0+' +lat_1='+lat1+' +lat_2='+lat2+' +x_0='+x0+' +y_0='+y0)
    return lcc
    

def tz_latlon(lat, lon, WHERETZ, JAN1, JUN1):
    """
    Timezone from latitude and longitude.
    :param lat: latitude [deg]
    :type lat: float
    :param lon: longitude [deg]
    :type lon: float
    :return: timezone
    :rtype: float
    """
    # Include module
    from tzwhere import tzwhere
    import pytz

    # get name of time zone using tzwhere, force to nearest tz
    #tz_name = WHERETZ.tzNameAt(lat, lon, forceTZ=True)
    tz_name = WHERETZ.tzNameAt(lat, lon)

    # check if coordinates are over international waters
    if not tz_name or tz_name in ('uninhabited', 'unknown'):
        # coordinates over international waters only depend on longitude
        return lon // 15.0, None
    else:
        tz_info = pytz.timezone(tz_name)  # get tzinfo

    # get the daylight savings time timedelta
    tz_date = JAN1  # standard time in northern hemisphere
    if tz_info.dst(tz_date):
        # if DST timedelta is not zero, then it must be southern hemisphere
        tz_date = JUN1  # a standard time in southern hemisphere
    tz_str = tz_info.localize(tz_date).strftime('%z')  # output timezone from ISO
    # convert ISO timezone string to float, including partial timezones
    return float(tz_str[:3]) + float(tz_str[3:]) / 60.0, tz_info

def scan_timezones(outfile, attr_in, *, loutf = False):
    # Include functions
    from CAMxtools.write.wrt_ioapi import wrt_ioapi
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools._cnvt._data2fin import _data2fin
    from tzwhere import tzwhere
    
    # Get nx, ny, dxy, and lcc
    nx = attr_in['NCOLS']
    ny = attr_in['NROWS']
    dxy = float(attr_in['XCELL']) # Simply assume XCELL = YCELL
    lcc = get_lcc(attr_in)

    # Global variables related to timezone
    # timezone lookup, force nearest tz for coords outside of polygons
    #global WHERETZ
    WHERETZ = tzwhere.tzwhere()
    # daylight savings time (DST) in northern hemisphere starts in March and ends
    # in November and the opposite in southern hemisphere
    #global JAN1
    JAN1 = datetime.datetime(2016, 1, 1)  # date with standard time in northern hemisphere
    #global JUN1
    JUN1 = datetime.datetime(2016, 6, 1)  # date with standard time in southern hemisphere
    
    # Calculating timezones over the domain
    print ("Calculating timezones over the domain")
    # If tz is set to auto, calculate a tshift array before the loop
    tzone_ji = np.zeros((ny,nx)).astype(int) # PST:-8 MST:-7 CST:-6 EST:-5
    for i in range ( nx ):
      for j in range ( ny ):
          lat, lon = proj_latlon_single(i,j,dxy,lcc)
          tzone_ji[j,i],tz_info = tz_latlon(lat, lon, WHERETZ, JAN1, JUN1) #tz_cell is based on LST not LDT.
    for itz in (np.unique(-tzone_ji)):
      print("time zone = {}".format(itz))

    # Data array preparation
    nspc = 1; nz = 1; nsteps = 1
    data2sav  = np.zeros((nspc,nsteps,nz,ny,nx))
    data2sav[0,0,0,:,:] = -tzone_ji
    tracernames = "TZONE".split()

    # Write output to a binary file
    fout = _data2fin(data2sav, tracernames, attr_in)
    l2uam = False # Force to output to IOAPI format
    if loutf: # Write output to netcdf
      if l2uam:
        wrt_uamiv(outfile, fout)
      else:
        wrt_ioapi(outfile, fout)
    else:
      return tracernames, data2sav

def main():
    # Include functions
    from CAMxtools.write.set_attr import set_attr

    ## define user inputs here ##
    outfile = str(sys.argv[1])
    infile  = str(sys.argv[2])

    # Check the first day file and read
    if not os.path.exists(infile):
      print ("{} does not exist!".format(infile))
      print ('Program exits from scan_timezones')
      exit()
    try:
      fin = uamiv(infile)
    except:
      try:
        fin = ncdf4.Dataset(infile)
      except:
        print ("Unrecognized file type")
        print ("infile = {}".format(infile))
        exit()

    # File attributes from an argument
    lfin0 = True
    fin0 = fin
    attr_in = set_attr(lfin0, fin0, {})
    loutf = True

    scan_timezones(outfile, attr_in, loutf = loutf)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
