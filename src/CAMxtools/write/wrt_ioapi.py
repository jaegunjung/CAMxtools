_ = r"""                                                                        
.. _dumper                                                                      
:mod:`dumper` -- CAMxtools dump module                                         
============================================                                    
                                                                                
... module:: wrt_ioapi
    :platform: Unix, Windows                                                    
    :synopsis: Write IOAPI file
    :details: Output IOAPI file using fin class which is defined
              in fclass in ../_cnvt/_data2fin.py
    :history:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>                              
"""                                               

__all__=['wrt_ioapi','find_unit',]
import netCDF4 as ncdf4
import datetime

def find_unit(spc_in,*,lemis=False,lcmaq=False):
    #unit for species
    CAMx_PMSPCs = "PNO3 PSO4 PNH4 POA SOA SOP PEC FPRM FCRS CPRM CCRS NA PCL PH2O".split() # microgram/m**3
    CMAQ_PMSPCs = "ASO4 ANH4 ANO3 AALK AXYL ATOL ABNZ ATRP AISO ASQT AORGC APOC APNCOM AEC AOTHR AFE AAL ASI ATI ACA AMG AK AMN ACORS ASOIL AH2O ANA ACL ASEACAT AISO3 AOLG".split() # microgram/m**3
    CMAQ_EMIS_PMSPCs = "PSO4 PNH4 PNO3 POC PNCOM PEC PMOTHR PFE PAL PSI PTI PCA PMG PK PMN PMC PH2O PNA PCL".split() # g/s
    CMAQ_SRFs = "SRFATKN SRFACC SRFCOR".split() # m**2/m**3
    CMAQ_NUMs = "NUMATKN NUMACC NUMCOR".split() # #/m**3

    if lemis:
      if lcmaq:
        unit = "moles/s"
        for spc in CMAQ_EMIS_PMSPCs:
          if spc_in.startswith(spc):
            unit = "g/s"
            break
      else:
        unit = "mole/hr"
        for spc in CAMx_PMSPCs:
          if spc_in.startswith(spc):
            unit = "g/hr"
            break
    else:
      unit = "ppmV"
      if lcmaq:
        for spc in CMAQ_PMSPCs:
          if spc_in.startswith(spc):
            unit = "microgram/m**3"
            break
        for spc in CMAQ_SRFs:
          if spc == spc_in:
            unit = "m**2/m**3"
            break
        for spc in CMAQ_NUMs:
          if spc == spc_in:
            unit = "#/m**3"
            break
      else:
        for spc in CAMx_PMSPCs:
          if spc_in.startswith(spc):
            unit = "microgram/m**3"
            break

    return unit

def wrt_ioapi(fout_path, fin, *, lsurf = False, lapp = False, lemis = False, lcmaq = False, ounits = None):
  """
  Create an IOAPI file using fin class.
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

  # Prepare file attributes
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

  # Open output file
  if lapp:
    try: fout = ncdf4.Dataset(fout_path,'a',format="NETCDF3_CLASSIC")
    except: exit("Output file does not exist while lapp is True")
  else:
    fout = ncdf4.Dataset(fout_path,'w',format="NETCDF3_CLASSIC")

  # Copy dimensions
  if not lapp:
    dimensions_keys = "TSTEP DATE-TIME LAY VAR".split() # These dimension keys must be in the file to be used by m3tools
    if lgrdded: # GRIDDED
      dimensions_keys.append("ROW"); dimensions_keys.append("COL")
    else:
      dimensions_keys.append("PERIM")
    for i in dimensions_keys:
      if i == 'VAR': size = novars
      elif i == 'TSTEP': size = None
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
      fout.createDimension(i,size=size)

  # Copy global attributes
  if not lapp:
    attribs = "IOAPI_VERSION EXEC_ID FTYPE CDATE CTIME WDATE WTIME SDATE STIME TSTEP NTHIK NCOLS NROWS NLAYS NVARS GDTYP P_ALP P_BET P_GAM XCENT YCENT XORIG YORIG XCELL YCELL VGTYP VGTOP VGLVLS GDNAM UPNAM VAR-LIST FILEDESC HISTORY".split() # These attributes must be in the file to be used by m3tools
    cdate = int(datetime.date.today().strftime("%Y%j"))
    ctime = int(datetime.datetime.now().strftime("%H%M%S"))
    for i in attribs:
      try: val = getattr(fin,i)
      except: val = ""
      if 'numpy.float32' in str(type(val)):
         val = val.item()
      if i == 'IOAPI_VERSION': val = '{:<80s}'.format("$Id: @(#) ioapi library version 3.1 $")
      if i == 'EXEC_ID': val = '{:<80s}'.format("wrt_ioapi in CAMxtools")
      if i == 'CDATE': val = cdate
      if i == 'CTIME': val = ctime
      if i == 'WDATE': val = cdate
      if i == 'WTIME': val = ctime
      if i == 'NTHIK': val = 1
      if i == 'NCOLS': val = nx
      if i == 'NROWS': val = ny
      if i == 'NLAYS':
         if lsurf:
           val = 1
         else:
           size = nz
      if i == 'NVARS': val = novars
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
      if i == 'GDNAM' or i == 'UPNAM':
         val = '{:<16s}'.format("unspecified")
      if i == 'FILEDESC': val = '{:<80s}'.format("unspecified")
      fout.setncattr(i,val)

  # Set start and end records of output file depending on lnew.
  start = 0
  tstep = getattr(fin,'TSTEP')
  if lapp: # Start = ("fin SDATE" - "fout SDATE") * 24000 / tstep
    start = int((datetime.datetime.strptime(str(getattr(fin,'SDATE')),"%Y%j")-datetime.datetime.strptime(str(getattr(fout,'SDATE')),"%Y%j")).days*240000/tstep)
  end = start + nsteps

  #copy variables
  extra = "TFLAG".split()
  for i in extra+varnames:
    if i == 'TFLAG':
      if not lapp:
        fout.createVariable(i,('int32'),(u'TSTEP', u'VAR', u'DATE-TIME'))
        fout.variables[i].setncattr("units","<YYYYDDD,HHMMSS>")
        fout.variables[i].setncattr("long_name","TFLAG")
        fout.variables[i].setncattr("var_desc",'{:<80s}'.format("Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS"))
      idate = getattr(fin,'SDATE')
      itime = getattr(fin,'STIME')
      for istep in range(start,end):
        fout.variables['TFLAG'][istep,:,0] = [idate for ivar in range(len(varnames))]
        fout.variables['TFLAG'][istep,:,1] = [itime for ivar in range(len(varnames))]
        itime += tstep
        if itime == 240000:
          itime = 0
          idate = int((datetime.datetime.strptime(str(idate),"%Y%j") + datetime.timedelta(days=1)).strftime("%Y%j"))
    else:
      if not lapp:
        if lgrdded: # GRIDDED
          if 'numpy.float' in str(type(fin.variables[var0][0,0,0,0])):
            dtype = 'float32'
          else:
            dtype = 'int32'
          fout.createVariable(i, (dtype), ('TSTEP', 'LAY', 'ROW', 'COL'))
        else:
          if 'numpy.float' in str(type(fin.variables[var0][0,0,0])):
            dtype = 'float32'
          else:
            dtype = 'int32'
          fout.createVariable(i, (dtype), ('TSTEP', 'LAY', 'PERIM'))
        #Find a relevant unit
        if lounit:
          s = varnames.index(i)
          unit = ounits[s]
        else:
          unit = find_unit(i,lemis=lemis,lcmaq=lcmaq)
        fout.variables[i].setncattr("long_name",'{:<16s}'.format(i))
        fout.variables[i].setncattr("units",'{:<16s}'.format(unit))
        fout.variables[i].setncattr("var_desc",'{:<80s}'.format("".join(["VARIABLE ",i])))
      if lgrdded: # GRIDDED
        fout.variables[i][start:end,:,:,:]=fin.variables[i][:,:,:,:]
      else:
        fout.variables[i][start:end,:,:]=fin.variables[i][:,:,:]

  # close files
  del fin
  fout.close()
  print ('*** SUCCESS writing netCDF file')
  return
