_ = r"""                                                                        
.. _dumper                                                                      
:mod:`dumper` -- CAMxtools dump module                                         
============================================                                    
                                                                                
... module:: calc_W126
    :platform: Unix, Windows                                                    
    :synopsis: Takes CAMx or CMAQ file header and tail names and calculates
               W126.
    :details: 1. Input files are CAMx or CMAQ binary output files.
              2. Output file is IOAPI formatted file with one time stamp.
    :history: - Corrected by adding checking lyyyyjjj to the first file to figure out file format (8/30/2017, jjung).
              - Based on python3 (8/7/2017, jjung)
              - Calculates 3-months sums and find maximum out of these 3-months sums (8/1/2017, jjung).
              - Output time stamp based on LST instead of UTC (5/3/2017, jjung).
              - Changes output format to IOAPI (5/1/2017, jjung).
              - Original version is prepared by lhuang (4/11/2017, lhuang).
              
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>                              
"""                                               

__all__=['calc_W126',]
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
import calendar
import numpy as np
import netCDF4 as ncdf4
import math
import datetime
import os

def calc_daily_W126 ( o3, nx, ny ) :
    W126 = np.zeros((ny,nx))
    for hr in range ( 12 ):
        for i in range ( nx ):
            for j in range ( ny ):
                    W126[j,i] = W126[j,i] + o3[hr,j,i] / ( 1 + 4403 * math.exp(-126 * o3[hr,j,i]))
    return W126

def calc_W126(outfile,fileh,filet,year,smo,emo,tz,*,lyyyyjjj=True,loutf=True,lverbose=True):
    # Include functions
    from CAMxtools.tzone.get_local_hrs import get_ozone_12hr
    from CAMxtools.write.wrt_ioapi import wrt_ioapi
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools._cnvt._data2fin import _data2fin

    # Check the first day file to figure out nx and ny
    jdate  = datetime.date(year,smo,1).strftime("%Y%j")
    gdate  = year*10000+smo*100+1
    if lyyyyjjj:
      infile = fileh+'.'+str(jdate)+'.'+filet
    else:
      infile = fileh+'.'+str(gdate)+'.'+filet
    if not os.path.exists(infile):
      print ("{} does not exist!".format(infile))
      print ('Program exits from calc_W126 at point 1')
      exit()
    try:
      fin = uamiv(infile)
      ftype = 'avg'
    except:
      try:
        fin = ncdf4.Dataset(infile)
        ftype = 'netcdf'
      except:
        print ("Unrecognized file type")
        print ("infile = {}".format(infile))
        exit()
    nx = len(fin.dimensions['COL'])
    ny = len(fin.dimensions['ROW'])
    
    # Loop to get monthly sums of W126
    nmo = emo - smo + 1
    monsum  = np.zeros((nmo,ny,nx))
    imo = 0
    mo = smo
    print ("Calculating monthly sum of W126")
    while imo < nmo:
      print ("Processing month :{}".format(calendar.month_name[mo]))
      jbdate  = int(datetime.date(year,mo,1).strftime("%Y%j"))
      edd = calendar.monthrange(year, mo)[1]
      jedate = int(datetime.date(year,mo,edd).strftime("%Y%j"))
      print ("jedate={}".format(str(jedate)))
      jdate = jbdate
      while (jdate <= jedate):
        jdatep1 = jdate + 1
        if lyyyyjjj:
            infile1 = fileh+'.'+str(jdate)+'.'+filet
            infile2 = fileh+'.'+str(jdatep1)+'.'+filet
        else:
            gdate = year*10000+mo*100+(jdate-jbdate+1)
            infile1 = fileh+'.'+str(gdate)+'.'+filet
            if (jdate < jedate):
                gdatep1 = gdate + 1
            else:
                gdatep1 = year*10000+(mo+1)*100+1
            infile2 = fileh+'.'+str(gdatep1)+'.'+filet
        
        for infile in [infile1, infile2]:
            if not os.path.exists(infile):
                print ("{} does not exist!".format(infile))
                print ('Program exits from calc_W126 at point 2')
                exit()
        
        if ftype == 'avg':
            fin1 = uamiv(infile1)
            fin2 = uamiv(infile2)
        else:
            fin1 = ncdf4.Dataset(infile1)
            fin2 = ncdf4.Dataset(infile2)
        
        # get ozone conc. from 8am-8pm local time
        if (lverbose): print ("  Processing Julian day :{}".format(str(jdate)))
        o3_48_utc = np.append(fin1.variables["O3"],fin2.variables["O3"],axis=0)
        o3 = get_ozone_12hr(o3_48_utc, tz, nx, ny)
    
        # calculate W126
        W126 = calc_daily_W126(o3 , nx, ny)
    
        # calculate monthly sum of W126
        for i in range ( nx ):
            for j in range ( ny ):
                monsum[imo,j,i] += W126[j,i]
        
        # increase day by one
        date = datetime.datetime.strptime(str(jdate),"%Y%j")
        date += datetime.timedelta(days=1)
        jdate = int(date.strftime("%Y%j"))
    
      # increase month by one
      imo += 1 # increase month counter
      mo += 1  # increase the order of month in the year
    
    # Loop to calculate three months sum
    nmo = emo - smo - 1 # remove two from the total number of months
    mon3sum  = np.zeros((nmo,ny,nx))
    imo = 0
    mo = smo
    print ("Calculating three monthly sum of W126")
    while imo < nmo:
      print ("Processing month :{}".format(calendar.month_name[mo]))
      for i in range ( nx ):
          for j in range ( ny ):
              mon3sum[imo,j,i] = monsum[imo,j,i] + monsum[imo+1,j,i] + monsum[imo+2,j,i]
      # increase month by one
      imo += 1 # increase month counter
      mo += 1  # increase the order of month in the year
    
    # Find the highest three monthly sum
    h1mon3sum  = np.zeros((ny,nx))
    imo = 0
    print ("Finding the highest three monthly sum of W126")
    while imo < nmo:
      for i in range ( nx ):
          for j in range ( ny ):
              if h1mon3sum[j,i] < mon3sum[imo,j,i]:
                  h1mon3sum[j,i] = mon3sum[imo,j,i]
      # increase month by one
      imo += 1 # increase month counter
      mo += 1  # increase the order of month in the year
    # write output to netcdf
    fin0 = fin # to copy file header
    nspc = 1; nz = 1; nsteps = 1
    data2sav  = np.zeros((nspc,nsteps,nz,ny,nx))
    data2sav[0,0,0,:,:] = h1mon3sum
    tracernames = "O3".split()
    beghr = 8

    #file attributes from an argument
    lfin0 = True
    attr_fed = {}
    attr_in = set_attr(lfin0, fin0, attr_fed , beghr = beghr)

    # Write output to a binary file
    fout = _data2fin(data2sav, tracernames, attr_in)
    l2uam = False # Force to output to IOAPI format
    if loutf:
      if l2uam:
        wrt_uamiv(outfile, fout)
      else:
        wrt_ioapi(outfile, fout)
    else:
      return tracernames, data2sav
