_ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: combine_timezones
    :platform: Unix, Windows
    :synopsis: Takes a csv file (inf_lst) which lists a set of IOAPI or 
               UAMIV files. Each file represents each time zone. Using
               IOAPI formatted tzfile which has timezone value at each
               grid cell, the set of files are chosen. The output file
               shows the same type of information as input files, but
               values corresponding time zones were chosen in the output
               file.
    :details:  The inf_lst has following format,
               # timezone,file
               6, /path1/file_cst.ncf
               7, /path2/file_mst.ncf
    :history:

... moduleauthor:: Jaegun Jung <jjung@ramboll.com>
"""

__all__=['combine_timezones',]
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
import numpy as np
import netCDF4 as ncdf4
import datetime
import os
import csv

def combine_timezones(outfile,inf_lst,tzfile,*,linfile=True,ltzfile=True,lverbose=True,attr_in0=None,nsteps=1,loutf=True):
    """
    Combine multiple data or files created based on each timezone into 
    relevant time at each grid cell by picking the relevant data/file
       - Inputs of concentration are binary files or data
       - Inputs of timezone are binary files or data
       - Output is a binary file
    Arguments:
       outfile  - IOAPI or UAMIV formatted output file
       inf_lst  - A list of files which have concentrations. If linfile
                 is set to False, it is data arrays instead of files
       tzfile   - IOAPI formatted binary file. If ltzfile is False, it
                 is an array.
       linfile  - False means inf_lst is a list of data array.
       ltzfile  - False means tzfile is a data array.
       lverbose - False means there is no screen outputs.
       attr_in0 - Input file attributes needed if linfile == False
       nsteps   - If three dimensionaly data feeded such as
                  data(nday,ny,nx), nsteps >= 2
       loutf    - Whether to return an output file (True) or not (False).
    """
    # Include functions
    from CAMxtools.tzone.get_local_hrs import get_ozone_12hr
    from CAMxtools.write.wrt_ioapi import wrt_ioapi
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools._cnvt._data2fin import _data2fin
         
    if ltzfile:
      print ("Reading the tzfile")
      print ("tzfile = {}".format(tzfile))
      ftz = ncdf4.Dataset(tzfile)
      tzone_ji = ftz.variables["TZONE"][0,0,:,:]
    else:
      tzone_ji = tzfile

    tzones = np.unique(tzone_ji)
    if lverbose:
      for itz in tzones:
        print("time zone = {}".format(itz))
    
    print ("Reading the input list text file")
    inf_tzs =[]
    if linfile:
      infiles =[]
      with open(inf_lst) as csvfile:
        lines = csv.reader(csvfile, delimiter=',', quotechar='|')
        for line in lines:
           if ''.join(line).strip() == '': continue # skip blank line
           if (line[0][0] == '/') or (line[0][0] == '#'): continue # skip the lines starting with /, #, or blank
           inf_tzs.append(int(line[0].split()[0]))
           infiles.append(line[1].split()[0])
    else:
      for itz in tzones:
        inf_tzs.append(int(itz))

    # Sanity check
    ntz = len(tzones)
    if linfile:
      no_infs = len(infiles)
    else:
      no_infs = len(inf_lst)
    if ntz != no_infs:
      print("A number of infiles does not match a number time zones in the time zone file")
      print("A number of infiles = {}".format(no_infs))
      print("A number of time zones from the time zone file = {}".format(ntz))
      exit()
    
    # Read input files or rename data
    if linfile: # Read input files
      print ("Reading 1st input file and diagnose")
      infile = infiles[0]
      fin = []
      if not os.path.exists(infile):
        print ("{} does not exist!".format(infile))
        print ('Program exits from combine_timezones')
        exit()
      try:
        fin0 = uamiv(infile)
        fin.append(fin0.variables["O3"][:,0,:,:])
        ftype = 'avg'
      except:
        try:
          fin0 = ncdf4.Dataset(infile)
          fin.append(fin0.variables["O3"][:,0,:,:])
          ftype = 'netcdf'
        except:
          print ("Unrecognized file type")
          print ("infile = {}".format(infile))
          exit()
      print ("Reading input files from 2nd one")
      for infile in infiles[1:]:
        print ("infile={}".format(infile))
        if ftype == 'avg':
            fin.append(uamiv(infile).variables["O3"][:,0,:,:])
        else:
            fin.append(ncdf4.Dataset(infile).variables["O3"][:,0,:,:])
      tz_fin_dict = dict(zip(inf_tzs,fin))
    else: # Rename data
      tz_fin_dict = dict(zip(inf_tzs,inf_lst)) # inf_lst are a list of data arrays
    
    print ("Finding the input file corresponding to time zone")
    if linfile:
      nx = len(ftz.dimensions['COL'])
      ny = len(ftz.dimensions['ROW'])
    else:
      nx = attr_in0['NCOLS']
      ny = attr_in0['NROWS']
    conc_lst = np.zeros((nsteps,ny,nx))
    for i in range ( nx ):
      for j in range ( ny ):
        if nsteps == 1:
          conc_lst[0,j,i] = tz_fin_dict[tzone_ji[j,i]][0,j,i]
        elif nsteps >= 2:
          conc_lst[:,j,i] = tz_fin_dict[tzone_ji[j,i]][:,j,i]
        else:
          print("Your no. of steps is {}".format(nsteps))
          exit("No of steps must be larger than 0")
    
    print ("Writing output file")

    # Data array preparation
    nspc = 1; nz = 1
    data2sav  = np.zeros((nspc,nsteps,nz,ny,nx))
    data2sav[0,:,0,:,:] = conc_lst
    tracernames = "O3".split()
    if linfile: # Create file attributes from fin0
      lfin0 = True
      attr_in = set_attr(lfin0, fin0, {})
    else:
      attr_in = attr_in0

    # Write output to a binary file
    fout = _data2fin(data2sav, tracernames, attr_in)
    l2uam = False # Force to output to IOAPI format
    if loutf: # Write output to netcdf
      if l2uam:
        wrt_uamiv(outfile, fout)
      else:
        wrt_ioapi(outfile, fout)
    else:
      return data2sav

def main():
    ## define user inputs here ##
    outfile = str(sys.argv[1])
    inf_lst = str(sys.argv[2])
    tzfile  = str(sys.argv[3])

    combine_timezones(outfile,inf_lst,tzfile)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
