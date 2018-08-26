__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: combine
    :platform: Unix, Windows
    :synopsis: Provide combine like functionality for CAMxtools
    :details: Takes either UAM or IOAPI as input file and outputs IOAPI files.
    :warning: 1. Don't use annual gigantic input files as this python script 
                 is very slow. Use the original combine program which comes 
                 with CMAQ.         
              2. In the species which starts with a number such as 0630101ALL. 
                 This type of variable name typically comes from CAMx DDM 
                 outputs. The current python script cannot handle this type
                 of variables.
    :history: 1. This version encapsulates (make a function) the contents in 
                 the main block to be referenced from other python program.
              2. Read readme_combine.txt for more details.
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['combine',]
from PseudoNetCDF.camxfiles.Memmaps import uamiv

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
import numpy as np
import netCDF4 as ncdf4
import csv
import re
import os
import datetime
import xarray as xr

def combine(outfile, comb, nifiles, infile, *, lsurf=True, lverbose=False, loutf = True, l2uam = False, lovarunits = False):
  """
  Linearly combine multiple variables in single or multiple input files
     - Input file can be UAM or IOAPI.
     - Output file format is IOAPI.
  Arguments:
     outfile - IOAPI or UAMIV formatted output file
     comb - species definition text file which has output species name, unit, and the equation of linear combination from species in input files. Follow the same syntax of CMAQ species definition file.
     nifiles - no. of input files
     infile - python list of multiple input files
     lsurf - if True, create only surface layer regardless of no. of layers from input file.
     lverbose - if True, print more messages to screen
     loutf - if True, create an output file
     l2uam - if True, output file format is UAM. Otherwise, IOAPI.
     lovarunits - if True, return output variable units, ovarunits.
  """
  # Include functions to call
  from CAMxtools.combine._exprparse import ExpressionEvaluator
  from CAMxtools.write.wrt_ioapi import wrt_ioapi
  from CAMxtools.write.wrt_uamiv import wrt_uamiv
  from CAMxtools.write.set_attr import set_attr
  from CAMxtools._cnvt._data2fin import _data2fin

  # Handle input files
  if lverbose: print ("no. of input files = {}".format(nifiles))
  fin = []
  for ifile in range(0,nifiles):
    try :
        fin.append(uamiv(infile[ifile]))
        if ifile == 0 : ftype = 'uam'
    except :
        try:
          fin.append(xr.open_dataset(infile[ifile]))
          if ifile == 0 : ftype = 'netcdf'
        except:
          print("Check whether your input file exists")
          print(infile[ifile])
          exit()
  
  # Used for variable names to recognize variables from equations
  NAME    = r'(?P<NAME>[a-zA-Z_][a-zA-Z_0-9\[\]]*)'

  # Outmost block begins with opening species definition file, comb.
  with open(comb) as csvfile:
      # Initialize tracernames which will be used for 'VAR-LIST' of output file at the end of this function.
      tracernames=[]
      ovarunits=[]

      # Scan the species defnition file to know no. of output variables.
      lines = csv.reader(csvfile, delimiter=',', quotechar='|')
      nspc = sum(1 for row in lines if not ((''.join(row).strip().replace(" ","") == '') or (row[0][0] == '/') or (row[0][0] == '#') or (row[0][0] == '!'))) # no. of output variables
      csvfile.seek(0) # rewind the species definition file

      # If an output file does not exist, get fin0 and nsteps
      lnew = False # flag for the new output
      if loutf:
         if not os.path.exists(outfile): lnew = True
      if ftype == 'uam':
         fin0 = fin[0]
         nsteps = len(fin[0].dimensions['TSTEP'])
      else:
         fin0 = ncdf4.Dataset(infile[0])
         nsteps = fin[0].dims['TSTEP']
      lfin0 = True
      attr_in = set_attr(lfin0, fin0, {})

      # Set data2sav
      lgrdded = False
      if attr_in['FTYPE'] == 1: lgrdded = True
      ny = attr_in['NROWS']
      nx = attr_in['NCOLS']
      nz = attr_in['NLAYS']
      if lsurf:
        nz = 1
      if lgrdded: # GRIDDED
        data2sav  = np.zeros((nspc,nsteps,nz,ny,nx))
      else: # BOUNDARY CONDITION
        ncells = 2*(nx+ny)+4
        data2sav  = np.zeros((nspc,nsteps,nz,ncells))

      # Main loop - process the species definition file line by line
      lines = csv.reader(csvfile, delimiter=',', quotechar='|')
      for line in lines:

          # Skip unnecessary lines
          if ''.join(line).strip().replace(" ","") == '': continue # skip blank line
          if ((line[0][0] == '/') or (line[0][0] == '#') or (line[0][0] == '!')): continue # skip the lines starting with /, #, or blank

          # Set ovar, ovarunit, and formula. Append tracernames.
          ovar = line[0].split()[0]
          ovarunit = line[1].split()[0]
          formula = line[2].split()[0]
          tracernames.append(ovar) 
          ovarunits.append(ovarunit) 

          # Find vars from formula, which are either an input file or ovar from previous lines, (i.e. [0])
          p = re.compile(NAME) # Declare pattern match that include alphabet and/or number and ends with "[ number ]"
          vars = p.findall(formula) # Find continuous blocks of p defined above. For example, if formula = 'NO_DD[1]+NO_DD[2]+NO2_DD[1]', vars = ['NO_DD[1]', 'NO_DD[2]', 'NO2_DD[1]']

          # Loop through individual var in vars
          for var in vars:
            # Delimit by [ or ] and store the first element to varname and the second element to fins. For example, varname[0] = NO_DD and fins[0] = 1. Do this for varname[1] and fins[1], and so on.
            varname = (re.findall(r"\w+",var)[0])
            findx = (int(re.findall(r"\w+",var)[1]))

            # If the var in formula is not defined, do sanity check fins and define. For example, NO_DD_1 = fin.variables[NO_DD][0:nstamps,:,:,:]
            var_findx = "".join([varname,"_",str(findx)])
            if not var_findx in locals():
              if findx < 0:
                print ('File index is negative! {}'.format(findx))
                exit()
              elif findx == 0: # if [0] in the var, use a variable that is already calculated.
                s = tracernames.index(varname)
                if lgrdded: # GRIDDED
                  exec("".join([var_findx, " = data2sav[s,:,0:nz,:,:]"]))
                else: # BOUNDARY CONDITION
                  exec("".join([var_findx, " = data2sav[s,:,0:nz,:]"]))
              elif findx <= nifiles:
                exec("".join([var_findx, " = fin[int(findx)-1].variables[varname][0:nsteps,0:nz,:,:]"]))
              else:
                print ('File index is larger than no. of input files! {}'.format(findx))
                print ('no. of input files {}'.format(nifiles))
                exit()

          # Change formular to an easy expression to deal with. For exmaple, NO_DD[1]+NO_DD[2]+NO2_DD[1] to NO_DD_1+NO_DD_2+NO2_DD_1
          for i in range(0,nifiles+1): # replace [] to _. (e.g. NO_DD[1] -> NO_DD_1)
            formula = formula.replace("".join(["[",str(i),"]"]),"".join(["_",str(i)]))

          # Construct dict_vars such as {'NO_DD_1':NO_DD_1, 'NO_DD_2':NO_DD_2, 'NO2_DD_1':NO2_DD_1}
          dict_vars={}
          vars = p.findall(formula) # return all the variables in formula
          for var in vars:
            exec("".join(["dict_vars['",var,"']=",var]))

          # Generate class e. To declare variables in dic_vars in the ExpressionEvaluator class. dic_vars is passed as an argument.
          e = ExpressionEvaluator(dict_vars)

          # Evaluate formula
          s = tracernames.index(ovar)
          if lgrdded: # GRIDDED
            data2sav[s,0:nsteps,0:nz,:,:] = e.parse(formula)
          else: # BOUNDARY CONDITION
            data2sav[s,0:nsteps,0:nz,:] = e.parse(formula)
  
      # Del files in fin lists
      if ftype == 'netcdf':
        for ifile in range(0,nifiles):
          fin[ifile].close()
      del fin

      # Return results
      if loutf: # if creating an output file
        fout = _data2fin(data2sav, tracernames, attr_in)
        lounit = True
        if l2uam:
          wrt_uamiv(outfile, fout, lsurf = lsurf, lapp = not lnew, ounits = ovarunits)
        else:
          wrt_ioapi(outfile, fout, lsurf = lsurf, lapp = not lnew, ounits = ovarunits)
      else:
        if lovarunits:
          return tracernames, data2sav, ovarunits
        else:
          return tracernames, data2sav

def main():
  # Check SURFACE_LAYER_ONLY environ variable
  lsurf = False
  try :
    surflag  = os.environ['SURFACE_LAYER_ONLY']
  except :
    surflag  = 'F'
  if (surflag.lower() == 't') or (surflag.lower() == 'y') :
     lsurf = True

  # Check OUT2UAM environ variable
  l2uam = False
  try :
    uamflag  = os.environ['OUT2UAM']
  except :
    uamflag  = 'F'
  if (uamflag == 'T') or (uamflag == 't') or (uamflag == 'Y') or (uamflag == 'y') :
     l2uam = True
  print ("Output only surface? {}".format(lsurf))
  outfile  = str(sys.argv[1])
  comb     = str(sys.argv[2])
  nifiles  = int(sys.argv[3]) # no. of input files
  infile = []
  for ifile in range(0,nifiles):
    infile.append(str(sys.argv[4+ifile]))
  
  loutf = True
  return combine(outfile, comb, nifiles, infile, lsurf=lsurf, loutf=loutf, l2uam=l2uam)

if __name__ == '__main__':
  # For internal use (no CAMxtools package installed), set the package path.
  try :
    package_path = os.environ['PACKAGE_PATH']
    sys.path.append(package_path)
  except :
    print ("PACKAGE_PATH environment variable is not set.")

  # Main
  main()
