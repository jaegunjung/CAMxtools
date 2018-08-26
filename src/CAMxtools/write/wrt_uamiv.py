_ = r"""                                                                        
.. _dumper                                                                      
:mod:`dumper` -- CAMxtools dump module                                         
============================================                                    
                                                                                
... module:: wrt_uamiv
    :platform: Unix, Windows                                                    
    :synopsis: Write UAMIV file
    :details: Output IOAPI file using fin class which is defined
              in fclass in ../_cnvt/_data2fin.py
    :history:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>                              
"""                                               

__all__=['wrt_uamiv',]

def wrt_uamiv(fout_path, fin, *, lsurf = False, lapp = False, lemis = False, lcmaq = False, ounits = None):
  """
  Create an UAMIV file using fin class.
  Arguments:
     fout_path - The IOAPI file name including path
     fin       - A class which include all the contents, data arrays, attributes, variable names.
                 The data array shape is (nspc,nsteps,nz,ny,nx).
     lsurf     - If True, output file has only 1 layer regardless of shape of data array.
     lapp      - If True, output file must already exist. Skip writing file dimensions and attributes.
                 It simply appends data arrays into the existing output file.
     lemis     - If True, it is emission file, and variable units have rates instead of concentration
     lcmaq     - If True, the file is for CMAQ, and varnames are checked against CMAQ PM species
                 so that PM unit is assigned, microgram/m**3 instead of moles/m**3.
     ounits    - A list which has units of variables
  """
  # Handling arguments
  lounit = True
  if ounits == None:
    lounit = False
    ounits = []

  # Include modules
  import datetime
  from PseudoNetCDF import PseudoNetCDFFile
  from PseudoNetCDF.pncgen import pncgen
  from PseudoNetCDF import PNC
  from PseudoNetCDF.sci_var import stack_files
  from CAMxtools.write.wrt_ioapi import find_unit
  from pathlib import Path
  import os

  #prepare file attributes
  novars = getattr(fin,'NVARS')
  varnames = getattr(fin,'VAR-LIST').split()
  var0 = varnames[0]
  nsteps = fin.variables[var0].shape[0]
  nz = getattr(fin,'NLAYS')
  ny = getattr(fin,'NROWS')
  nx = getattr(fin,'NCOLS')
  lgrdded = False
  if getattr(fin,'FTYPE') == 1: lgrdded = True
  if lgrdded: # GRIDDED
    assert len(fin.variables[var0].shape) == 4, "len(fin.variables[var0].shape) MUST be 4"
  else: # BOUNDARY CONDITION
    assert len(fin.variables[var0].shape) == 3, "len(fin.variables[var0].shape) MUST be 3"
    ncells = fin.variables[var0].shape[2]
    assert ncells == 2*(nx+ny)+4, "ncells MUST be 2*(nx+ny)+4"

  #open output file
  if not lgrdded: #BOUNDARY
    ftype = 'lateral_boundary'
  else:
    ftype = 'uamiv'
  newf = PseudoNetCDFFile()

  # Set tstep
  tstep = getattr(fin,'TSTEP')

  #copy dimensions
  dimensions_keys = "TSTEP DATE-TIME LAY VAR".split() # These dimension keys must be in the file to be used by m3tools
  if lgrdded: # GRIDDED
    dimensions_keys.append("ROW"); dimensions_keys.append("COL")
  else:
    dimensions_keys.append("PERIM")
  for i in dimensions_keys:
    if i == 'VAR': size = novars
    elif i == 'TSTEP': size = nsteps
    elif i == 'DATE-TIME': size = 2
    elif i == 'LAY':
       if lsurf:
         size = 1
       else:
         size = nz
    elif i == 'VAR': size = novars
    elif i == 'ROW': size = ny
    elif i == 'COL': size = nx
    else : size = ncells # i == 'PERIM'
    newf.createDimension(i,size)

  #copy global attributes
  attribs = "XORIG YORIG XCELL YCELL PLON PLAT TLAT1 TLAT2 IUTM ISTAG CPROJ GDTYP XCENT YCENT P_ALP P_BET P_GAM NLAYS NROWS NCOLS NVARS VAR-LIST NAME NOTE ITZON FTYPE VGTYP VGTOP VGLVLS GDNAM UPNAM FILEDESC SDATE STIME TSTEP Conventions history".split() # These attributes must be in the file for pncgen to uamiv
  cdate = int(datetime.date.today().strftime("%Y%j"))
  ctime = int(datetime.datetime.now().strftime("%H%M%S"))
  for i in attribs:
    try: val = getattr(fin,i)
    except: val = ""
    if 'numpy.float32' in str(type(val)):
       val = val.item()
    if i == 'PLON': val = getattr(fin,'XCENT')
    if i == 'PLAT': val = getattr(fin,'YCENT')
    if i == 'TLAT1': val = getattr(fin,'P_ALP')
    if i == 'TLAT2': val = getattr(fin,'P_GAM')
    if i == 'IUTM':
      val = 0
      if getattr(fin,'GDTYP') == 5: val = 1
    if i == 'ISTAG': val = 0
    if i == 'CPROJ':
      if getattr(fin,'GDTYP') == 1: #LATLON
        val = 0
      elif getattr(fin,'GDTYP') == 5: #UTM
        val = 1
      elif getattr(fin,'GDTYP') == 2: #LCP
        val = 2
      elif getattr(fin,'GDTYP') == 6: #PSP
        val = 4
      elif getattr(fin,'GDTYP') == 7: #Equatorial Mercator
        val = 5
      else:
        print("GDTYP = {}".format(GDTYP))
        exit("Not relevant projection")
    if i == 'NLAYS':
       if lsurf:
         val = 1
       else:
         size = nz
    if i == 'NROWS': val = ny
    if i == 'NCOLS': val = nx
    if i == 'NVARS': val = novars
    if i == 'NSTEPS': val = nsteps
    if i == 'NAME':
      name_str = 'AVERAGE'
      if lemis: name_str = 'EMISSIONS'
      if not lgrdded: name_str = 'BOUNDARY'
      val = '{:<10s}'.format(name_str)
    if i == 'NOTE': val = '{:<60s}'.format("wrt_uamiv in CAMxtools")
    if i == 'ITZON': val = 0
    if i == 'VGTYP':
       if lsurf or nz == 1:
         val = -9999
       else:
         size = 2 # VGSGPN3 non-h sigma-p
    if i == 'VGTOP':
       if lsurf or nz == 1:
         val = -9.999E36
       else:
         size = getattr(fin,'VGTOP')
    if i == 'GDNAM' or i == 'UPNAM' or i == 'FILEDESC':
       val = '{:<16s}'.format("CAMx")
    if i == 'Conventions': val = "CF-1.6"
    if i == 'history': val = '{:<250s}'.format("unspecified")
    setattr(newf,i,val)

  #copy variables
  extra = "TFLAG".split()
  for i in extra+varnames:
    if i == 'TFLAG':
      newf.createVariable(i,('int32'),(u'TSTEP', u'VAR', u'DATE-TIME'))
      newf.variables[i].setncattr("units","<YYYYDDD,HHMMSS>")
      newf.variables[i].setncattr("long_name","TFLAG")
      newf.variables[i].setncattr("var_desc",'{:<80s}'.format("Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS"))
      idate = getattr(fin,'SDATE')
      itime = getattr(fin,'STIME')
      for istep in range(nsteps):
        newf.variables['TFLAG'][istep,:,0] = [idate for ivar in range(len(varnames))]
        newf.variables['TFLAG'][istep,:,1] = [itime for ivar in range(len(varnames))]
        itime += tstep
        if itime == 240000:
          itime = 0
          idate = int((datetime.datetime.strptime(str(idate),"%Y%j") + datetime.timedelta(days=1)).strftime("%Y%j"))
    else:
      if lgrdded: # GRIDDED
        var_type = str(type(fin.variables[var0][0,0,0,0]))
      else:
        var_type = str(type(fin.variables[var0][0,0,0]))
      if 'numpy.float' in var_type:
        dtype = 'float32'
      else:
        dtype = 'int32'
      if lgrdded: # GRIDDED
        newf.createVariable(i, (dtype), ('TSTEP', 'LAY', 'ROW', 'COL'))
      else:
        newf.createVariable(i, (dtype), ('TSTEP', 'LAY', 'PERIM'))
      #Find a relevant unit
      if lounit:
        s = varnames.index(i)
        unit = ounits[s]
      else:
        unit = find_unit(i,lemis=lemis,lcmaq=lcmaq)
      newf.variables[i].setncattr("long_name",'{:<16s}'.format(i))
      newf.variables[i].setncattr("units",'{:<16s}'.format(unit))
      newf.variables[i].setncattr("var_desc",'{:<80s}'.format("".join(["VARIABLE ",i])))
      newf.variables[i]=fin.variables[i]

  # Create UAMIV file
  if lapp: # Save the original fout_path to fout_old_path
    fout_old_path = fout_path + ".old"
    try:
      os.rename(fout_path,fout_old_path) # Move the original output to fout_old_path
      #pncargs = '--format=' + ftype + ',mode="r+"'
      #pncargs = '--format=' + ftype
    except:
      print("Output file is {}".format(fout_path))
      exit("Output file does not exist while lapp is True")
    pncargs = '--format=' + ftype
    oldfile = PNC(pncargs, fout_old_path).ifiles[0]
    fout_tmp_path = fout_path + ".tmp"
    pncgen(newf, fout_tmp_path, format = ftype)
    tmpfile = PNC(pncargs, fout_tmp_path).ifiles[0]
    newf_app = stack_files([oldfile, tmpfile], 'TSTEP')
    pncgen(newf_app, fout_path, format = ftype)
    if Path(fout_old_path).exists(): os.remove(fout_old_path)
    if Path(fout_tmp_path).exists(): os.remove(fout_tmp_path)
  else:
    pncgen(newf, fout_path, format = ftype)

  # close files
  if 'newfile' in locals(): del newfile
  del fin
  if 'newf' in locals(): newf.close()
  print ('*** SUCCESS writing UAMIV file')
  return
