__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: get_local_hours
    :platform: Unix, Windows
    :synopsis: This is a collection of functions that extracts
          subset of hours from two consecutive days after
          converting to local standard time.
    :details:  Functions are, 
           1. get_local_24hr
           2. get_ozone_12hr
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['get_local_24hr','get_ozone_12hr','get_mda1s','get_mda8s','get_mda8s_a0','get_mda1','get_mda8','get_davg','get_dbav3']

def get_local_24hr (tracernames, indata, tz, ny, nx) :
    """
    Get local 24 hour data from two consecutive files based on UTC
    Arguments:
       tracernames - species names to extract
       indata - input data based on UTC for two consecutive days
       tz - array for local time zone (EST = 5, CST = 6, ...)
       ny - no. of rows in input file
       nx - no. of cols in input file
    """
    # Include functions
    import numpy as np

    ntracers = len(tracernames)
    trc4d_hrly = np.zeros ((ntracers,24,ny,nx))
    for j in range(ny):
      for i in range(nx):
        trc4d_hrly[:,:,j,i] = indata[:,tz[j,i]:tz[j,i]+24,0,j,i]
    del indata

    return trc4d_hrly

def get_ozone_12hr ( o3_48_utc, tz, nx, ny ) :
     # Include functions
     import numpy as np

     o3 = np.zeros((12,ny,nx))
     for j in range(ny):
       for i in range(nx):
         o3[:,j,i] = o3_48_utc[tz[j,i]+8:tz[j,i]+20,0,j,i]

     del o3_48_utc

     return o3

def get_mda1s ( conc_hrs_utc, tz, nd, nx, ny ) :
     # Include functions
     import numpy as np
     from PseudoNetCDF.userfuncs import daymax

     mda1s = np.zeros((nd,ny,nx))
     for j in range(ny):
       for i in range(nx):
         mda1s[:,j,i] = daymax ( conc_hrs_utc[tz[j,i]:,j,i], axis=0)

     del conc_hrs_utc

     return mda1s

def get_mda8s_a0 ( conc_hrs_utc, tz, nd, nx, ny ) :
     """
     Count only 17 Avg8hrs starting from 7 AM LST following a new definition of MDA8 from EPA.
     It returns multiple mda8s.
     """
     # Include functions
     import numpy as np
     from CAMxtools._psdncdf.userfuncs_a0 import mda8 as calc_mda8_a0

     mda8 = np.zeros((nd,ny,nx))
     for j in range(ny):
       for i in range(nx):
         mda8[:,j,i] = calc_mda8_a0 ( conc_hrs_utc[tz[j,i]:,j,i], axis = 0)

     del conc_hrs_utc

     return mda8

def get_mda8s ( conc_hrs_utc, tz, nd, nx, ny ) :
     """
     Old method by counting 24 Avg8hrs
     It returns multiple mda8s.
     """
     # Include functions
     import numpy as np
     from PseudoNetCDF.userfuncs import mda8 as calc_mda8

     mda8 = np.zeros((nd,ny,nx))
     for j in range(ny):
       for i in range(nx):
         mda8[:,j,i] = calc_mda8 ( conc_hrs_utc[tz[j,i]:,j,i], axis = 0)

     del conc_hrs_utc

     return mda8

def get_mda1 ( concs_48_utc, tz, nx, ny ) :
     # Include functions
     import numpy as np

     nspc = concs_48_utc.shape[0]
     mda1 = np.zeros((nspc,ny,nx))
     ind_hr = np.zeros((nspc,ny,nx)).astype(int)
     ind_hr_1spc = np.zeros((ny,nx)).astype(int)
     for j in range(ny):
       for i in range(nx):
         for ispc in range(nspc):
           concs_24_lst_1spc_1cell = concs_48_utc[ispc,tz[j,i]:tz[j,i]+24,j,i]
           ind_hr_1spc[j,i] = np.argsort(concs_24_lst_1spc_1cell, axis = 0)[concs_24_lst_1spc_1cell.shape[0]-1,...]
           ind_hr[ispc,j,i] = ind_hr_1spc[j,i]
           mda1[ispc,j,i] = concs_48_utc[ispc,tz[j,i]+ind_hr[ispc,j,i],j,i]

     del concs_48_utc

     return mda1, ind_hr

def get_mda8 ( concs_48_utc, tz, nx, ny, *, lnew_mda8 = False ) :
     # Include functions
     import numpy as np

     nspc = concs_48_utc.shape[0]
     mda8 = np.zeros((nspc,ny,nx))
     ind_hr = np.zeros((nspc,ny,nx)).astype(int)
     ind_hr_1spc = np.zeros((ny,nx)).astype(int)
     for j in range(ny):
       for i in range(nx):
         for ispc in range(nspc):
           concs_8hravg_utc_1spc_1cell = np.apply_along_axis(np.convolve, axis = 0, arr = concs_48_utc[ispc,tz[j,i]:tz[j,i]+31,j,i], v = [1/8.]*8, mode = 'valid') # 31 = 24 + 7
           nhrs2skip = 0
           if lnew_mda8: nhrs2skip = 7
           ind_hr_1spc[j,i] = np.argsort(concs_8hravg_utc_1spc_1cell[nhrs2skip:], axis = 0)[concs_8hravg_utc_1spc_1cell.shape[0]-nhrs2skip-1,...]
           ind_hr[ispc,j,i] = ind_hr_1spc[j,i]+nhrs2skip #based on LST
           mda8[ispc,j,i] = concs_8hravg_utc_1spc_1cell[ind_hr[ispc,j,i]]

     del concs_48_utc

     return mda8, ind_hr

def get_davg ( concs_48_utc, tz, nx, ny ) :
     # Include functions
     import numpy as np

     nspc = concs_48_utc.shape[0]
     davg = np.zeros((nspc,ny,nx))
     for j in range(ny):
       for i in range(nx):
         davg[:,j,i] = np.mean ( concs_48_utc[:,tz[j,i]:tz[j,i]+24,j,i], axis = 1)

     del concs_48_utc

     return davg

def get_dbav3 ( concs_48_utc, tz, nx, ny ) :
     # Include functions
     import numpy as np

     nspc = concs_48_utc.shape[0]
     dbav3 = np.zeros((8,nspc,ny,nx))
     for j in range(ny):
       for i in range(nx):
         for b in range(8):
           dbav3[b,:,j,i] = np.average ( concs_48_utc[:,tz[j,i]+b*3:tz[j,i]+b*3+3,j,i], axis = 1)

     del concs_48_utc

     return dbav3
