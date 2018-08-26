__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: hxmdaz_naaqs
    :platform: Unix, Windows
    :synopsis: It reads CAMx or CMAQ output files and prepares
               grand average for the period specified.
    :details:  The sequential processes are as follows,
                1. PROCESSING COMBINE: creates indata
                2. PROCESSING MDA8O3 or MDA1O3
                3. FIND 1st or 4th highest or average of MDAZO3
                4. CREATE CMAQ or CAMx output file with one time stamp:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['hxmdaz_naaqs',]
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

def hxmdaz_naaqs(jdbeg,jdend,outfile,csvfile,comb,attr_in,nifiles,infile_h,infile_t,*,avg_hr=8,rank=4,lyyyyjjj=True,tzone=None,l2uam=False,lnew_mda8=False,tmpoutf=None):
    # Include functions
    from CAMxtools.combine.combine import combine
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr
    from CAMxtools.tzone.scan_timezones import scan_timezones
    from CAMxtools.tzone.get_local_hrs import get_mda1s
    from CAMxtools.tzone.get_local_hrs import get_mda8
    from CAMxtools._cnvt._data2fin import _data2fin
    from CAMxtools.write.wrt_ioapi import wrt_ioapi
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.write.wrt_uamiv import wrt_uamiv
    from CAMxtools.naaqs.wrt_csv_for_naaqs import wrt_csv_for_naaqs
    import numpy as np
    import gc

    # Check arguments
    l1tzone = False
    if tzone != None: l1tzone = True # If users specify a specific time zone, MATS MDA8 O3 is calculated based on the time zone.

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

    # Get attribute
    nx  = attr_in['NCOLS']
    ny  = attr_in['NROWS']

    # Set a variable name 
    spec = "O3"

    # Daily loop
    jdate = jdbeg
    while (jdate <= jdend):
      print ("Processing {}".format(jdate))
      gdate = int(datetime.datetime.strptime(str(jdate),"%Y%j").strftime("%Y%m%d"))
      infile = []
      for ifile in range(0,nifiles):
        if lyyyyjjj:
          infile.append(infile_h[ifile]+'.'+str(jdate)+'.'+infile_t[ifile])
        else:
          infile.append(infile_h[ifile]+'.'+str(gdate)+'.'+infile_t[ifile])
      print ("  1. PROCESSING COMBINE")
      tracernames, indata2, ovarunits = combine(None,comb,nifiles,infile,lverbose=False,loutf=False,lovarunits=True)
      if (jdate > jdbeg):
        # One time execution - Check 'O3' is the first output species and index output trcnames which excludes O3 from tracernames
        if jdate == jdbeg + 1:
          s = tracernames.index(spec)
          if s != 0: exit('O3 must be the first species in combine spec_def')
          ntracers = len(tracernames)
          trcnames = tracernames[1:]
        print ("  2. PROCESSING METRICS FOR PREVIOUS DAY")
        concs_48_utc = np.append(indata1[:,:,0,:,:],indata2[:,:,0,:,:],axis=1)
        conc_o3_48_utc = concs_48_utc[np.newaxis,s,:,:,:] #concs_48_utc[nspc,nt,ny,nx], conc_o3_48_utc[1,nt,ny,nx]
        if avg_hr == 1:
          dmdaz_1day, dind_hr_1day = get_mda1(conc_o3_48_utc, tzone_ji, nx, ny) # dmda1_1day[1,ny,nx], dind_hr_1day[1,ny,nx]
        elif avg_hr == 8:
          dmdaz_1day, dind_hr_1day = get_mda8(conc_o3_48_utc, tzone_ji, nx, ny, lnew_mda8 = lnew_mda8) # dmda1_1day[1,ny,nx], dind_hr_1day[1,ny,nx]
        else:
          exit("avg_hr must be either 1 or 8.")
        # Set trcs_dmdaz_1day
        trcs_dmdaz_1day = np.zeros((ntracers-1,ny,nx))
        for i in range(nx):
          for j in range(ny):
            if avg_hr == 1:
              trcs_dmdaz_1day[:,j,i] = concs_48_utc[1:,np.squeeze(dind_hr_1day)[j,i],j,i]
            elif avg_hr == 8:
              trcs_daz_1day_1cell = np.apply_along_axis(np.convolve, axis = 1, arr = concs_48_utc[1:,tzone_ji[j,i]:tzone_ji[j,i]+31,j,i], v = [1/8.]*8, mode = 'valid') #trcs_daz_1day_1cell[nspc-1,24]
              trcs_dmdaz_1day[:,j,i] = trcs_daz_1day_1cell[:,np.squeeze(dind_hr_1day)[j,i]]
            else:
              exit("avg_hr must be either 1 or 8.")
          gc.collect()
        if jdate == jdbeg + 1:
          dmdazs = dmdaz_1day #dmdazs[nd,ny,nx]
          dind_hrs = dind_hr_1day #dind_hrs[nd,ny,nx]
          trcs_dmdazs = np.array([trcs_dmdaz_1day]) #trcs_dmdazs[nd,nspc-1,ny,nx]
        if jdate > jdbeg + 1:
          dmdazs = np.append(dmdazs,dmdaz_1day,axis=0)
          dind_hrs = np.append(dind_hrs,dind_hr_1day,axis=0)
          trcs_dmdazs = np.append(trcs_dmdazs,np.array([trcs_dmdaz_1day]),axis=0)
      indata1 = indata2
      jdate = int((datetime.datetime.strptime(str(jdate),"%Y%j") + datetime.timedelta(days=1)).strftime("%Y%j"))
      gc.collect()
    del indata2
    del indata1
    gc.collect()

    # Prepare MDAZ if asked
    lout_mdaz = False
    if not tmpoutf == None: lout_mdaz = True
    if lout_mdaz:
      nd = jdend - jdbeg
      data2sav = np.zeros((ntracers,nd,1,ny,nx))
      data2sav[0,:,0,:,:] = dmdazs
      data2sav[1:,:,0,:,:] = np.einsum('jikl->ijkl',trcs_dmdazs)
      # Write a binary file for MDAZ
      attr_in['TSTEP']=240000
      fout = _data2fin(data2sav, tracernames, attr_in)
      if l2uam:
        wrt_uamiv(tmpoutf, fout, lsurf = True, ounits = ovarunits)
      else:
        wrt_ioapi(tmpoutf, fout, lsurf = True, ounits = ovarunits)
      gc.collect()

    # Find X highest
    if rank == 1:
      rank_name = "FIRST"
    elif rank == 4:
      rank_name = "FOURTH"
    elif rank == 0:
      rank_name = "AVERAGE"
    else:
      exit("rank must be either 1 or 4 or 0.")
    if rank == 0:
      print ("  3. PROCESSING AVERAGE")
    else:
      print ("  3. PROCESSING {} HIGHEST".format(rank_name))
    data2sav = np.zeros((ntracers,1,1,ny,nx))
    nd = jdend - jdbeg
    if (nd < rank):
      print('Your no. of days is less than the rank you specified.')
      print('No. of days from input files = {}'.format(nd))
      print('Rank you select is = {}'.format(rank))
      exit('Either reduce your rank or increase no. of days from input files')
    if rank == 0:
      data2sav[0,0,0,:,:] = np.mean(dmdazs,axis=0) #data2sav[nspc,nt,nz,ny,nx], dmdazs[nd,ny,nx]
      data2sav[1:,0,0,:,:] = np.mean(trcs_dmdazs,axis=0) #trcs_dmdazs[nd,nspc-1,ny,nx]
    else:
      ind = np.argsort(dmdazs,axis=0)[nd-rank,...] #dmdazs[nd,ny,nx], ind[ny,nx]
      gr = np.ogrid[0:dmdazs.shape[0],0:dmdazs.shape[1],0:dmdazs.shape[2]]
      gr[0] = ind
      data2sav[0,0,0,:,:] = dmdazs[gr] #data2sav[nspc,nt,nz,ny,nx]
      for i in range(nx):
        for j in range(ny):
          data2sav[1:,0,0,j,i] = trcs_dmdazs[ind[j,i],:,j,i]
    # Write a csv file
    if rank == 0:
      jdays = np.zeros((ny,nx))+attr_in['SDATE'] # Average of MDAZ, set the beginning date
      hours = np.zeros((ny,nx)) # Average of MDAZ, set zeros.
    else:
      jdays = attr_in['SDATE']+ind[:,:]
      hours = np.zeros((ny,nx))+dind_hrs[gr][0,:,:] #dind_hrs[gr].shape = (1,ny,nx)
    data2csv = np.zeros((ntracers,ny,nx))
    data2csv = data2sav[:,0,0,:,:]
    wrt_csv_for_naaqs( csvfile, tracernames, data2csv, jdays, hours)
    gc.collect()
    # Write a binary file for HXMDAZ
    fout = _data2fin(data2sav, tracernames, attr_in)
    if l2uam:
      wrt_uamiv(outfile, fout, lsurf = True, ounits = ovarunits)
    else:
      wrt_ioapi(outfile, fout, lsurf = True, ounits = ovarunits)
    gc.collect()

def main():
    # Include functions to call
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr

    # Check TZONE environment variable
    try :
      tz = int(os.environ['TIMEZONE'])
    except :
      tz = None
    # Check OUT2UAM environment variable
    l2uam = False
    try :
      uamflag  = os.environ['OUT2UAM']
    except :
      uamflag  = 'F'
    if (uamflag == 'T') or (uamflag == 't') or (uamflag == 'Y') or (uamflag == 'y') :
       l2uam = True
    # Check MDAZ_OUTF environment variable
    try :
      tmpoutf  = os.environ['MDAZ_OUTF']
    except :
      tmpoutf = None
    if not tmpoutf == None:
       print("  MDAZ_OUTF will be created as {}".format(tmpoutf))

    # Check arguments
    outfile  = str(sys.argv[1])
    csvfile  = str(sys.argv[2])
    avg_hr   = int(sys.argv[3])
    rank     = int(sys.argv[4])
    jdbeg    = int(sys.argv[5])
    jdend    = int(sys.argv[6])
    lyyyyjjj = False
    if str(sys.argv[7]).lower() == "true" :
      lyyyyjjj = True # CAMx file name has yyyyjjj
    comb     = str(sys.argv[8])
    nifiles  = int(sys.argv[9])
    infile_h = []
    infile_t = []
    # Build input file header and tail lists
    nargs2drop = 10 # No. of arguments up to "infile1 header" (inclusive)
    for ifile in range(0,nifiles):
      infile_h.append(str(sys.argv[nargs2drop+(2*ifile)]))
      infile_t.append(str(sys.argv[nargs2drop+1+(2*ifile)]))

    # Check NEW_MDA8_METHOD environment variable
    lnew_mda8 = False
    try :
      new_mda8_flag  = os.environ['NEW_MDA8_METHOD']
    except :
      new_mda8_flag  = 'F'
    if (new_mda8_flag == 'T') or (new_mda8_flag == 't') or (new_mda8_flag == 'Y') or (new_mda8_flag == 'y') :
       lnew_mda8 = True
    else:
      if avg_hr == 8:
        print("  WARNING: You ARE USING OLD MDA8 METHOD!\n 24 Running 8 hour averages will be considered to calculate MDA8 O3 instead of 17 averages")

    # Check the first day file and read
    jdate  = str(jdbeg)
    gdate  = int(datetime.datetime.strptime(str(jdate),"%Y%j").strftime("%Y%m%d"))
    if lyyyyjjj:
      infile0 = infile_h[0] + '.' + str(jdate) + '.' + infile_t[0]
    else:
      infile0 = infile_h[0] + '.' + str(gdate) + '.' + infile_t[0]
    if not os.path.exists(infile0):
      print ("{} does not exist!".format(infile0))
      print ('Program exits from hxmdaz_naaqs')
      exit()
    try:
      fin0 = uamiv(infile0)
    except:
      try:
        fin0 = ncdf4.Dataset(infile0)
      except:
        print ("Unrecognized file type")
        print ("infile1 = {}".format(infile0))
        exit()
    
    # Set file attributes
    lfin0 = True
    attr_fed = {}
    attr_in = set_attr(lfin0,fin0,attr_fed)

    # Do the main process
    hxmdaz_naaqs(jdbeg,jdend,outfile,csvfile,comb,attr_in,nifiles,infile_h,infile_t,avg_hr=avg_hr,rank=rank,lyyyyjjj=lyyyyjjj,tzone=tz,l2uam=l2uam,lnew_mda8=lnew_mda8,tmpoutf=tmpoutf)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
