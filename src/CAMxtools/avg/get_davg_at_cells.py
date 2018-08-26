__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: get_davg_at_cells
    :platform: Unix, Windows
    :synopsis: It calculates daily average at Class I or II areas
               based on local time stamp.
    :details:  The sequential processes are as follows,
                1. Read two consecutive CAMx or CMAQ daily files
                   whose timestep is an hour.
                2. Do sanity check whether the two files have
                   same no. of columns and rows.
                3. Applying timeshift (tzone), construct 24 hour
                   data based on local time from 48 hour data
                   provided based on UTC.
                4. Calculate daily average.
                5. Write to csv file at Class I or II areas using
                   wrt_csv_for_vis.
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['get_davg_at_cells',]
from PseudoNetCDF.camxfiles.Memmaps import uamiv
import numpy as np
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
import netCDF4 as ncdf4
import os

def get_davg_at_cells(tracernames,indata,outfile,xref,yyyyjjj,*,tzone=None):
    # Include functions
    from CAMxtools.vis.wrt_csv_for_vis import wrt_csv_for_vis
    from CAMxtools.tzone.get_local_hrs import get_local_24hr
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools._cnvt._data2fin import _data2fin

    # Begin the script
    ny = indata.shape[3]
    nx = indata.shape[4]
    ntracers = len(tracernames)
    
    # read the first variable and set the dimension
    trc4d_hrly   = get_local_24hr(tracernames, indata, tzone, ny, nx)
    
    # define the variables
    trc3d_davg_ntrcs = np.zeros((ntracers,ny,nx))
    
    # calculate the daily averages
    trc3d_davg_ntrcs[:,:,:] = np.average(trc4d_hrly[:,:,:,:],axis=1)
    
    #   write the csv file
    wrt_csv_for_vis (tracernames, trc3d_davg_ntrcs, xref, yyyyjjj, outfile)

def main():
    # Include functions to call
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools.tzone.scan_timezones import scan_timezones
    import datetime

    # CHECK USER INPUTS
    # Check TZONE environment variables
    try :
      tzone = int(os.environ['TIMEZONE'])
    except:
      tzone = None
    # Check TZONE environment variables
    try :
      tzfile = str(os.environ['TZFILE'])
    except:
      tzfile = None
    # If a user specified a time zone, do not use tzfile.
    if not tzone == None:
      tzfile = None # If a user specified a time zone, do not use tzfile.
    # Arguments
    outfile = str(sys.argv[1])
    infile1 = str(sys.argv[2])
    infile2 = str(sys.argv[3])
    xref    = str(sys.argv[4])

    # Begin the script
    try:
      fin1 = uamiv(infile1)
      fin2 = uamiv(infile2)
    except:
      try:
        fin1 = ncdf4.Dataset(infile1)
        fin2 = ncdf4.Dataset(infile2)
      except:
        print ("Unrecognized file type")
        print ("infile1 = {}".format(infile1))
        print ("infile2 = {}".format(infile2))
        exit()
    tracernames = list(fin1.variables.keys())
    tracernames.remove('TFLAG')

    nx = len(fin1.dimensions['COL'])
    ny = len(fin1.dimensions['ROW'])
    nx2 = len(fin2.dimensions['COL'])
    ny2 = len(fin2.dimensions['ROW'])
    if (nx != nx2) or (ny != ny2):
      print ("infile1 and infile2 have different no. of columns or rows")
      print ("infile1 x = {}, y = {}".format(nx,ny))
      print ("infile2 x = {}, y = {}".format(nx2,ny2))
      exit()

    # Prepare tzone 2d data array
    lfin0 = True
    attr_fed = {}
    attr_in = set_attr(lfin0,fin1,attr_fed)
    l1tzone = False
    if tzone != None:
      l1tzone = True # If users specify a specific time zone, vis is calculated based on the time zone.
    if not tzfile == None:
      print ("Reading the tzfile")
      print ("tzfile = {}".format(tzfile))
      ftz = ncdf4.Dataset(tzfile)
      tzone_ji = ftz.variables["TZONE"][0,0,:,:].astype(np.int)
    else:
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

    # Set indata
    ntracers = len(tracernames)
    nt = len(fin1.dimensions['TSTEP'])
    nz = attr_in['NLAYS']
    indata1 = np.zeros((ntracers,nt,nz,ny,nx))
    indata2 = np.zeros((ntracers,nt,nz,ny,nx))
    for varname in tracernames:
      s = tracernames.index(varname)
      indata1[s,:,:,:,:] = fin1.variables[varname][:,:,:,:]
      indata2[s,:,:,:,:] = fin2.variables[varname][:,:,:,:]
    indata = np.append(indata1,indata2,axis=1)

    # Main process
    yyyyjjj = getattr(fin1,'SDATE')
    get_davg_at_cells(tracernames,indata,outfile,xref,yyyyjjj,tzone=tzone_ji)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
