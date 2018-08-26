__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: regrid
    :platform: Unix, Windows
    :synopsis: It reads CAMx or CMAQ output files and regrids
               to the nested domain.
    :details:  The sequential processes are as follows,
                1. PROCESSING COMBINE: creates indata2
                2. PROCESSING REGRID:
                3. CREATE OUTPUT FILE:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['distance','closest_ij','find_ij','regrid',]
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

def distance(lonA,latA,lonB,latB):
    """
    Calculates a distance between point A to point B and
    returns it.
    Arguments:
       lonA - longitude at the point A
       latA - latitude at the point A
       lonB - longitude at the point B
       latB - latitude at the point B
    """
    # include modules
    from math import sin, cos, sqrt, atan2, radians

    # Radius of earth in km
    try:
      R = float(os.environ['IOAPI_ISPH'])/1000.
    except:
      R = 6370.0

    # Radius of earth in km
    lat1 = radians(latA); lon1 = radians(lonA)
    lat2 = radians(latB); lon2 = radians(lonB)

    # Radius of earth in km
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    d = R * c

    return d

def closest_ij(icol_min,icol_max,jcol_min,jcol_max,lons_in,lats_in,lon_out,lat_out):
    """
    Find a closest (i,j) of an input domain for (lon_out,lat_out).
    returns it.
    Arguments:
       icol_min - a minimum column index for a column loop of the input domain
       icol_max - a maximum column index for a column loop of the input domain
       jcol_min - a minimum row index for a column loop of the input domain
       jcol_max - a maximum row index for a column loop of the input domain
       ny_in - no. of rows in the input domain
       lons_in - a 2D array of longitudes for the input domain
       lats_in - a 2D array of latitudes for the input domain
       lon_out - longitude at one point whose closest lons_in and lats_in need to be found
       lat_out - latitude at one point whose closest lons_in and lats_in nee
d to be found
    """
    min_d = 9.99E36
    iloc_1cell = 99999; jloc_1cell = 99999
    for j in range(jcol_min,jcol_max+1):
        for i in range(icol_min,icol_max+1):
            d = distance(lons_in[j,i],lats_in[j,i],lon_out,lat_out)
            if min_d > d:
                min_d = d
                iloc_1cell = i; jloc_1cell = j
    
    return iloc_1cell, jloc_1cell

def find_ij(dxy_in, x0, y0, nx_in, ny_in, nx_out, ny_out, xpos, ypos, *, lno_edges = True,lnearest=False,lats_in=None,lons_in=None,lats_out=None,lons_out=None):
    """
    Receives a parent domain definition and returns icell and jcell
    indexes corresponding to a nest domain.
    Arguments:
       dxy_in - delta x or y of the parent domain
       x0     - x value at the SW corner of the parent domain
       y0     - y value at the SW corner of the parent domain
       nx_in  - no. of columns of the parent domain
       ny_in  - no. of rows of the parent domain
       nx_out - no. of columns of the nest domain
       ny_out - no. of rows of the nest domain
       xpos   - the parent domain's x value at the center of 
                the nest domain's grid cells. It has the nest
                domain's 2d array shape.
       ypos   - the parent domain's y value at the center of 
                the nest domain's grid cells. It has the nest
                domain's 2d array shape.
       lno_edges - Exclude buffer cells from the parent domain
                   if True.
    """
    # Include functions
    import numpy as np

    # Assign the range of grid cell index of parent domain
    icol_min = 1; icol_max = nx_in-2
    jcol_min = 1; jcol_max = ny_in-2
    if not lno_edges:
      icol_min = 0; icol_max = nx_in-1
      jcol_min = 0; jcol_max = ny_in-1

    # Main loop
    iloc = np.zeros((ny_out,nx_out),dtype=int); jloc = np.zeros((ny_out,nx_out),dtype=int)
    for j in range ( ny_out ) :
        for i in range ( nx_out ) :
            iloc[j,i] = int((xpos[j,i] - x0)/dxy_in)
            jloc[j,i] = int((ypos[j,i] - y0)/dxy_in)
            if iloc[j,i] < icol_min or iloc[j,i] > icol_max or jloc[j,i] < jcol_min or jloc[j,i] > jcol_max:
              if lnearest:
                iloc[j,i], jloc[j,i] = closest_ij(icol_min,icol_max,jcol_min,jcol_max,lons_in,lats_in,lons_out[j,i],lats_out[j,i])
              else:
                print('Nested grid extends outside parent grid')
                print('West-East parent range        : {}, {}'.format(icol_min,icol_max))
                print('West-East nest cell location  : {}'.format(iloc[j,i]))
                print('North-South parent range      : {}, {}'.format(jcol_min,jcol_max))
                print('North-South nest cell location: {}'.format(jloc[j,i]))
                exit('Make the nested grid smaller')
            print('i,j,iloc,jloc= {}, {}, {}, {}'.format(i,j,iloc[j,i],jloc[j,i]))
    return iloc, jloc

def regrid(outfile,project,utmzon,plon,plat,tlat1,tlat2,xorg,yorg,dxy_out,nx_out,ny_out,comb,attr_in,nstep,nifiles,infile,*,lno_edges=True,l2uam=False,lnearest=False,ijfile=None):
    # Include functions
    from CAMxtools.combine.combine import combine
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools.tzone.scan_timezones import get_lcc
    from CAMxtools.regrid.projection import ll2latlon
    from CAMxtools.regrid.projection import lcp2latlon
    from CAMxtools.regrid.projection import ll2lcp
    from CAMxtools._cnvt._data2fin import _data2fin
    from CAMxtools.write.wrt_ioapi import wrt_ioapi
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    import numpy as np
    import csv

    # 1. COMBINE
    print ("  1. PROCESSING COMBINE")
    tracernames, indata, ovarunits = combine(None,comb,nifiles,infile,lsurf=False,lverbose=False,loutf=False,lovarunits=True)

    # 2. REGRID
    print ("  2. PROCESSING REGRID")

    # 2.1 Set attributes of output file
    if l2uam:
      NA_val = 0.
    else:
      NA_val = -9.999E36
    attr_out = {}
    for key, value in attr_in.items():
      if key == 'XCELL': attr_out[key] = dxy_out
      elif key == 'YCELL': attr_out[key] = dxy_out
      elif key == 'XORIG': attr_out[key] = xorg
      elif key == 'YORIG': attr_out[key] = yorg
      elif key == 'NCOLS': attr_out[key] = nx_out
      elif key == 'NROWS': attr_out[key] = ny_out
      elif key == 'GDTYP':
        if project == 'LATLON': attr_out[key] = 1
        elif project == 'LAMBERT': attr_out[key] = 2
        else: exit()
      elif key == 'P_ALP':
        if project == 'LATLON': attr_out[key] = NA_val
        elif project == 'LAMBERT': attr_out[key] = tlat1
        else: exit()
      elif key == 'P_BET':
        if project == 'LATLON': attr_out[key] = NA_val
        elif project == 'LAMBERT': attr_out[key] = tlat2
        else: exit()
      elif key == 'P_GAM':
        if project == 'LATLON': attr_out[key] = NA_val
        elif project == 'LAMBERT': attr_out[key] = plon
        else: exit()
      elif key == 'XCENT':
        if project == 'LATLON': attr_out[key] = NA_val
        elif project == 'LAMBERT': attr_out[key] = plon
        else: exit()
      elif key == 'YCENT':
        if project == 'LATLON': attr_out[key] = NA_val
        elif project == 'LAMBERT': attr_out[key] = plat
        else: exit()
      else:
        attr_out[key] = value

    if ijfile == None:
      # 2.2a.1 Get lat and lon at center of grid cells for output
      if attr_out['GDTYP'] == 1:
        lats_out, lons_out = ll2latlon (xorg, yorg, dxy_out, nx_out, ny_out, lno_edges = False)
      elif attr_out['GDTYP'] == 2:
        lcc_out = get_lcc(attr_out)
        lats_out, lons_out = lcp2latlon (lcc_out, dxy_out, nx_out, ny_out, lno_edges = False)
      else:
        print("Your GDTYP is {}".format(attr_in['GDTYP']))
        exit("This program currently supports LATLON (GDTYP = 1) and LCP (GDTYP = 2) only.")
      # 2.2a.2 Get x and y values of input file corresponding to each output grid cell
      if attr_in['GDTYP'] == 1:
        xpos = lons_out; ypos = lats_out
      elif attr_in['GDTYP'] == 2:
        x0_in = attr_in['XORIG']
        y0_in = attr_in['YORIG']
        lcc_in = get_lcc(attr_in)
        xpos, ypos = ll2lcp (lons_out, lats_out, lcc_in, ny_out, nx_out, x0_in, y0_in)
      else:
        print("Your GDTYP is {}".format(attr_in['GDTYP']))
        exit("This program currently supports LATLON (GDTYP = 1) and LCP (GDTYP = 2) only.")
      # 2.2a.3 Get grid index of input file corresponding to each output grid cell
      dxy_in = attr_in['XCELL']
      x0     = attr_in['XORIG']
      y0     = attr_in['YORIG']
      nx_in  = attr_in['NCOLS']
      ny_in  = attr_in['NROWS']
      if lnearest:
        if attr_in['GDTYP'] == 1:
          lats_in, lons_in = ll2latlon (x0, y0, dxy_in, nx_in, ny_in, lno_edges = False)
        elif attr_in['GDTYP'] == 2:
          lats_in, lons_in = lcp2latlon (lcc_in, dxy_in, nx_in, ny_in, lno_edges = False)
        else:
          print("Your GDTYP is {}".format(attr_in['GDTYP']))
          exit("This program currently supports LATLON (GDTYP = 1) and LCP (GDTYP = 2) only.")
        iloc, jloc = find_ij(dxy_in, x0, y0, nx_in, ny_in, nx_out, ny_out, xpos, ypos, lno_edges=True,lnearest=lnearest,lats_in=lats_in,lons_in=lons_in,lats_out=lats_out,lons_out=lons_out)
      else:
        iloc, jloc = find_ij(dxy_in, x0, y0, nx_in, ny_in, nx_out, ny_out, xpos, ypos, lno_edges=True)
    else:
      # 2.2b Read index mapping from ijfile
      iloc = np.zeros((ny_out,nx_out),dtype=int); jloc = np.zeros((ny_out,nx_out),dtype=int)
      with open(ijfile) as csvfile:
        lines = csv.reader(csvfile, delimiter=',', quotechar='|') # skip the first header line
        for line in lines:
          if line[0][0] == 'i': continue # skip the first header line
          i = int(line[0].split()[0])
          j = int(line[1].split()[0])
          iloc[j,i] = int(line[2].split()[0])
          jloc[j,i] = int(line[3].split()[0])

    # 2.5 Mapping and prepare data to save
    nspc = len(tracernames)
    nz = attr_out['NLAYS']
    data2sav = np.zeros((nspc,nstep,nz,ny_out,nx_out))
    for j in range (ny_out):
      for i in range (nx_out):
        data2sav[:,:,:,j,i] = indata[:,:,:,jloc[j,i],iloc[j,i]]
    del indata

    # 3. WRITE A BINARY FILE
    fout = _data2fin(data2sav, tracernames, attr_out)
    if l2uam:
      wrt_uamiv(outfile, fout, lsurf = False, ounits = ovarunits)
    else:
      wrt_ioapi(outfile, fout, lsurf = False, ounits = ovarunits)

def main():
    # Include functions to call
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr

    # Check INCLUDE_BUFFER_CELLS environment variable
    lno_edges = True
    try:
      bfflag = os.environ['INCLUDE_BUFFER_CELLS']
    except:
      bfflag = 'F'
    if (bfflag.lower() == 't') or (bfflag.lower() == 'y'):
      lno_edges = False
    # Check OUT2UAM environment variable
    l2uam = False
    try :
      uamflag  = os.environ['OUT2UAM']
    except :
      uamflag  = 'F'
    if (uamflag == 'T') or (uamflag == 't') or (uamflag == 'Y') or (uamflag == 'y') :
       l2uam = True
    # Check NEAREST environment variable
    lnearest = False
    try:
      nearflag = os.environ['NEAREST_CELL']
    except:
      nearflag = 'F'
    if (nearflag.lower() == 't') or (nearflag.lower() == 'y'):
      lnearest = True
    # Check IJFILE environment variable
    ijfile = None
    try:
      ijfile = os.environ['IJFILE']
    except:
      ijfile = None

    # Check arguments
    outfile = str(sys.argv[1])
    project = str(sys.argv[2]) # Output projection
    proj_list = ['LAMBERT', 'LATLON']
    if not project in proj_list:
      exit('Available output projections are {}.'.format(proj_list))
    utmzon  = int(sys.argv[3])    # UTM zone, MUST be 0 for now.
    plon    = float(sys.argv[4])  # Center/Pole lat
    plat    = float(sys.argv[5])  # Center/Pole lon
    tlat1   = float(sys.argv[6])  # True lat1
    tlat2   = float(sys.argv[7])  # True lat2
    xorg    = float(sys.argv[8])  # XORG
    yorg    = float(sys.argv[9])  # YORG
    dxy     = float(sys.argv[10]) # DX or DY
    nx      = int(sys.argv[11])   # NX
    ny      = int(sys.argv[12])   # NY
    comb    = str(sys.argv[13])
    nifiles = int(sys.argv[14])
    infile = []
    # Build an infile list
    nargs2drop = 15 # No. of arguments up to "infile1 header" (inclusive)
    for ifile in range(0,nifiles):
      infile.append(str(sys.argv[nargs2drop+ifile]))

    # Check the first file and read
    infile0 = infile[0]
    if not os.path.exists(infile0):
      print ("{} does not exist!".format(infile0))
      print ('Program exits from regrid')
      exit()
    try:
      fin0 = uamiv(infile0)
      nstep = len(fin0.dimensions['TSTEP'])
    except:
      try:
        fin0 = ncdf4.Dataset(infile0)
        nstep = fin0.dims['TSTEP']
      except:
        print ("Unrecognized file type")
        print ("infile1 = {}".format(infile0))
        exit()

    # If NEAREST_CELL is on and INCLUDE_BUFFER_CELLS is on, check 
    spc0 = getattr(fin0,'VAR-LIST').split()[0]
    val0 = fin0.variables[spc0][0,0,0,0]
    if lnearest:
      if lno_edges:
        if val0 != 0.0: print("WARNING: The value at the SW corner cell is not zero while you want to exclude buffer cells.")
      else:
        assert val0 != 0.0, "The value at the SW corner cell is zero while you want to include buffer cells. Program exits."

    # Set file attributes
    lfin0 = True
    attr_fed = {}
    attr_in = set_attr(lfin0,fin0,attr_fed)

    # Do the main process
    regrid(outfile,project,utmzon,plon,plat,tlat1,tlat2,xorg,yorg,dxy,nx,ny,comb,attr_in,nstep,nifiles,infile,lno_edges=lno_edges,l2uam=l2uam,lnearest=lnearest,ijfile=ijfile)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
