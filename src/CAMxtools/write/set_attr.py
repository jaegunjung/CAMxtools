_ = r"""                                                                        
.. _dumper                                                                      
:mod:`dumper` -- CAMxtools dump module                                         
============================================                                    
                                                                                
... module:: set_attr
    :platform: Unix, Windows                                                    
    :synopsis: Set global attributes of IOAPI file
    :details: The global attributes set from this function are,
                FTYPE, SDATE, STIME, TSTEP, NCOLS*,
                NROWS*, NLAYS**, GDTYP, P_ALP, P_BET, P_GAM, XCENT,
                YCENT, XORIG, YORIG, XCELL, YCELL,
                VGTYP, VGTOP, VGLVLS
              *NCOLS, NROW, NLAYS need to be set in case
              output array size is to be decided
              **NLAYS is set as it can be used related to lsurf
              flag from caller before calling wrt_ioapi.

              Following attributes are not set as they are either
              trivial or real data need to be known.
                CDATE, CTIME, WDATE, WTIME, NCOLS,
                NROWS, NVARS, GDNAM, UPNAM, VAR-LIST, 
                FILEDESC, HISTORY
              If these attributes are for boundary condition
              (FTYPE == 2), NCOLS and NROWS must be specified
              as the real data shape does not provide this
              information.
    :warning:
    :history:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>                              
"""                                               

__all__=['set_attr',]

def set_attr(lfin0,fin0,attr_fed,*,beghr = 0):
  """
  Set global attributes for an IOAPI file to write
     attr_in - returning dictionary which has global attributes
  Arguments:
     lfin0 - logical variable whether to use fin0 to get attributes
     fin0 - file already opened whose global attributes to be copied if lfin0 == True
     attr_fed - dictionary to use for the global attributes if lfin0 == False
     beghr - STIME global attribute can be overridden if provided
  """
  attr_in = {}
  if not lfin0:
    attr_in['FTYPE'] = attr_fed['FTYPE'] # 1 = GRIDDED, 2 = BOUNDARY
    attr_in['SDATE'] = attr_fed['SDATE'] #YYYYJJJ
    attr_in['STIME'] = attr_fed['STIME']
    try:
      attr_in['TSTEP'] = attr_fed['TSTEP']
    except:
      attr_in['TSTEP'] = 10000
    attr_in['NCOLS'] = attr_fed['NCOLS']
    attr_in['NROWS'] = attr_fed['NROWS']
    attr_in['NLAYS'] = attr_fed['NLAYS']
    if attr_fed['FTYPE'] == 2:
      attr_in['NCOLS'] = attr_fed['NCOLS']
      attr_in['NROWS'] = attr_fed['NROWS']
    attr_in['GDTYP'] = attr_fed['GDTYP'] # 1 = LATGRD3, 2 = LAMGRD3, 3 = MERGRD3, 4 = STEGRD3, 5 = UTMGRD3, 6 = POLGRD3, 7 = EQMGRD3, 8 = TRMGRD3, 9 = ALBGRD3, 10 = LEQGRD3 (See /models/CMAQ/lib64/ioapi-3.1/ioapi/fixed_src/PARMS3.EXT for more details)
    attr_in['P_ALP'] = attr_fed['P_ALP']
    attr_in['P_BET'] = attr_fed['P_BET']
    attr_in['P_GAM'] = attr_fed['P_GAM']
    attr_in['XCENT'] = attr_fed['XCENT']
    attr_in['YCENT'] = attr_fed['YCENT']
    attr_in['XORIG'] = attr_fed['XORIG']*1000.
    attr_in['YORIG'] = attr_fed['YORIG']*1000.
    attr_in['XCELL'] = attr_fed['XCELL']*1000.
    attr_in['YCELL'] = attr_fed['YCELL']*1000.
    attr_in['VGTYP'] = attr_fed['VGTYP'] # 1 = VGSGPH3, 2 = VGSGPN3 (non-h sigma-P), 3 = VGSIGZ3, 4 = VGPRES3, 5 = VGZVAL3, 6 = VGHVAL3, 7 = VGWRFEM, 8 = VGWRFNM (See /models/CMAQ/lib64/ioapi-3.1/ioapi/fixed_src/PARMS3.EXT for more details)
    attr_in['VGTOP'] = attr_fed['VGTOP']*100 # Pa
    attr_in['VGLVLS'] = attr_fed['VGLVLS']
  else:
    attribs = "FTYPE SDATE STIME TSTEP NCOLS NROWS NLAYS GDTYP P_ALP P_BET P_GAM XCENT YCENT XORIG YORIG XCELL YCELL VGTYP VGTOP VGLVLS".split() # These attributes must be in the file to be used by m3tools
    for i in attribs:
      if i == 'STIME':
        if beghr == 0:
          attr_in[i] = getattr(fin0,i)
        else:
          attr_in[i] = beghr * 10000
      else:
        attr_in[i] = getattr(fin0,i)
    del fin0
  return attr_in
