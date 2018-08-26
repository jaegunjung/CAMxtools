__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: get_vis_at_cells
    :platform: Unix, Windows
    :synopsis: It is a main script that processes visibility analysis from
               CAMx or CMAQ output files to (1) delta deci-view (DDV) of max,
               8th, W20%, B20%, and min over one year at each class 1 or 2
               area and (2) no. of days whose DDVs are larger than 0.5 or 1.0 
               DDV at each class 1 or 2.
    :details: In the main function, followings are sequentially processed,
                1. Sanity check arguments fed through a shell script
                2. Build input file header and tail lists used for combine.py.
                   The first file is indexed [1], and the second file is indexed
                   [2], and so on in the species definition file in combine.py.
                   (Note this is not for index of day. Date in the input files
                   will be fed in the day_loop_for_vis function.)
                3. "day_loop_for_vis" function runs for a year looping by day.
                4. "maxvis_temporal" function processes and reports DDV of max,
                   8th, W20%, B20%, and min at each class 1 or 2
                5. "count_temporal" function processes and reports no. of days
                   whose DDV is larger than 0.5 or 1 at each class 1 or 2.
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['get_vis_at_cells',]
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

def get_vis_at_cells(jdbeg,jdend,outf_ddv,outf_cnt,xref,comb,nifiles,infile_h,infile_t,class_lst,*,lyyyyjjj=True,tzone=None,attr_in=None):

    # Include functions to call
    import numpy as np
    from CAMxtools.vis.maxvis_temporal import maxvis_temporal
    from CAMxtools.vis.count_temporal import count_temporal
    from CAMxtools.vis.day_loop_for_vis import day_loop_for_vis
    from CAMxtools.tzone.scan_timezones import scan_timezones

    # Check arguments
    l1tzone = False
    if tzone != None:
      l1tzone = True # If users specify a specific time zone, vis is calculated based on the time zone.
      if attr_in == None:
        exit("When tzone is None, i.e auto-tzone, attr_in must not be None")

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

    # Run for a year looping by day
    allday_visf  = "allday_visf.csv"
    day_loop_for_vis(jdbeg,jdend,allday_visf,xref,comb,nifiles,infile_h,infile_t,class_lst,lyyyyjjj=True,tzone=tzone_ji)

    # Process and report DDV of max, 8th, W20%, B20%, and min at each class 1 or 2
    maxvis_temporal(allday_visf,outf_ddv)

    # Process and report no. of days whose DDV is larger than 0.5 or 1 at each class 1 or 2
    count_temporal(allday_visf,outf_cnt)

    # Remove the csv file that has annual visiblity results
    os.remove(allday_visf)

def list_args_exit(nifiles,nargs2drop,*,lyyyyjjj=True):
    print('')
    print('='*70)
    print('                           ARGUMENTS PROVIDED')
    print('='*70)
    print('Outfile = {}'.format(str(sys.argv[1])))
    print('Combine species definition file = {}'.format(str(sys.argv[2])))
    print('Class I and II Cross reference file = {}'.format(str(sys.argv[3])))
    print('Class I and II list file to map vis data = {}'.format(str(sys.argv[4])))
    print('Begining Julian Date = {}'.format(int(sys.argv[5])))
    print('Ending Julian Date = {}'.format(int(sys.argv[6])))
    print('Use Julian Date in input file = {}'.format(lyyyyjjj))
    print('No. of input files = {}'.format(int(sys.argv[8])))
    print('Input file for the first day')
    for ifile in range(0,nifiles):
      if lyyyyjjj:
        print('  input file{} = {}'.format(ifile+1,str(sys.argv[nargs2drop+(2*ifile)])+'.'+str(sys.argv[6])+'.'+str(sys.argv[nargs2drop+1+(2*ifile)])))
      else:
        gdbeg = int(datetime.datetime.strptime(str(sys.argv[6]),"%Y%j").strftime("%Y%m%d"))
        print('  input file{} = {}'.format(ifile+1,str(sys.argv[nargs2drop+(2*ifile)])+'.'+str(gdbeg)+'.'+str(sys.argv[nargs2drop+1+(2*ifile)])))
    print('='*70)
    exit()

def main():
    # Include functions to call
    import netCDF4 as ncdf4
    from PseudoNetCDF.camxfiles.Memmaps import uamiv
    from CAMxtools.write.set_attr import set_attr

    # CHECK USER INPUTS
    nargs2drop = 9 # No. of arguments upto "infile1 header" (inclusive)
    # Check TZONE and SURFACE_LAYER environment variables
    try :
      tz = int(os.environ['TIMEZONE'])
      if (tz < 0) or (tz > 12):
        print('CURRENT VERSION ONLY SUPPORT 0 <= time zone <= 12.') 
        list_args_exit(int((len(sys.argv)-nargs2drop)/2),nargs2drop,lyyyyjjj=lyyyyjjj)
    except :
      tz = None 
    lsurf = False
    try :
      surflag = os.environ['SURFACE_LAYER_ONLY']
    except :
      surflag = 'F'
    if (surflag == 'T') or (surflag == 't') or (surflag == 'Y') or (surflag == 'y') :
       lsurf = True
    print ("Output only surface? {}".format(lsurf))
    # Check arguments
    outname  = str(sys.argv[1])
    outf_ddv = outname + '.ddv.csv'
    outf_cnt = outname + '.cnt.csv'
    if Path(outf_ddv).exists(): os.remove(outf_ddv)
    if Path(outf_cnt).exists(): os.remove(outf_cnt)
    # Read lyyyyjjj first as it is an argument for list_args_exit
    lyyyyjjj = False
    if str(sys.argv[7]).lower() == "true" :
      lyyyyjjj = True # CAMx file name has yyyyjjj
    comb     = str(sys.argv[2])
    if not (Path(comb).exists()):
      print('COMBINE SPECIES DEFINITION FILE DOES NOT EXIST!')
      list_args_exit(int((len(sys.argv)-nargs2drop)/2),nargs2drop,lyyyyjjj=lyyyyjjj)
    xref = str(sys.argv[3])
    if not (Path(xref).exists()):
      print('CROSS REFERENCE FILE DOES NOT EXIST!')
      list_args_exit(int((len(sys.argv)-nargs2drop)/2),nargs2drop,lyyyyjjj=lyyyyjjj)
    class_lst = str(sys.argv[4])
    if not (Path(xref).exists()):
      print('CLASS12 FILE MAP TO IMPROVE SITE FOR VISIBILITY DOES NOT EXIST!')
      list_args_exit(int((len(sys.argv)-nargs2drop)/2),nargs2drop,lyyyyjjj=lyyyyjjj)
    jdbeg = int(sys.argv[5])
    try:
      datetime.datetime.strptime(str(jdbeg),"%Y%j")
    except:
      print('BEGINNING JULIAN DATE IS INVALID!')
      list_args_exit(int((len(sys.argv)-nargs2drop)/2),nargs2drop,lyyyyjjj=lyyyyjjj)
    jdend = int(sys.argv[6])
    try:
      datetime.datetime.strptime(str(jdend),"%Y%j")
    except:
      print('ENDING JULIAN DATE IS INVALID!')
      list_args_exit(int((len(sys.argv)-nargs2drop)/2),nargs2drop,lyyyyjjj=lyyyyjjj)
    if (jdbeg > jdend):
      print('BEGINNING JULIAN IS LARGER THAN ENDING!')
      list_args_exit(int((len(sys.argv)-nargs2drop)/2),nargs2drop,lyyyyjjj=lyyyyjjj)
    nifiles  = int(sys.argv[8]) # no. of input files
    infile_h = []
    infile_t = []
    if nifiles < (len(sys.argv)-nargs2drop)/2:
      print('NO. OF INPUT FILES IS SET LESS THAN NO. OF INPUT FILE ARGUMENTS PROVIDED!')
      list_args_exit(int((len(sys.argv)-nargs2drop)/2),nargs2drop,lyyyyjjj=lyyyyjjj)
    if nifiles > (len(sys.argv)-nargs2drop)/2:
      print('NO. OF INPUT FILES IS SET LARGER THAN NO. OF INPUT FILE ARGUMENTS PROVIDED!')
      list_args_exit(int((len(sys.argv)-nargs2drop)/2),nargs2drop,lyyyyjjj=lyyyyjjj)

    # Build input file header and tail lists
    for ifile in range(0,nifiles):
      infile_h.append(str(sys.argv[nargs2drop+(2*ifile)]))
      infile_t.append(str(sys.argv[nargs2drop+1+(2*ifile)]))

    # Check the first day file and read
    jdate  = str(jdbeg)
    gdate  = int(datetime.datetime.strptime(str(jdate),"%Y%j").strftime("%Y%m%d")) 
    if lyyyyjjj:
      infile0 = infile_h[0] + '.' + str(jdate) + '.' + infile_t[0]
    else:
      infile0 = infile_h[0] + '.' + str(gdate) + '.' + infile_t[0]
    if not os.path.exists(infile0):
      print ("{} does not exist!".format(infile0))
      print ('Program exits from run_calc_W126')
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

    # File attributes from an argument
    lfin0 = True
    attr_fed = {}
    attr_in = set_attr(lfin0,fin0,attr_fed)

    # Get annual visibility analysis results at Class I and II areas
    get_vis_at_cells(jdbeg,jdend,outf_ddv,outf_cnt,xref,comb,nifiles,infile_h,infile_t,class_lst,lyyyyjjj=lyyyyjjj,tzone=tz,attr_in=attr_in)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path. 
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
