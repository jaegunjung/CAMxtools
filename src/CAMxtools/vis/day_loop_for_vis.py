__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: day_loop_for_vis
    :platform: Unix, Windows
    :synopsis: It is a main body of get_vis_at_cells. Taking CAMx or
               CMAQ output files, calculate delta deci-view(DDV) every
               day at each grid cell such as class 1 or 2. This result
               is accumulated to allday vis file csv file, allday_visf.
    :details:  The sequential processes are as follows,
                1. PROCESSING COMBINE: creates indata2 which has
                       species for visility calculation such as NO2, _NH4_2SO4,
                       NH4NO3, EC, OM, Soil, CM, and SeaSalt
                2. PROCESSING DAILY AVERAGE FOR PREVIOUS DAY: creates
                       tmp_davg_out.csv, tmp_davg_outf.
                3. PROCESSING DECI-VIEW FOR PREVIOUS DAY: creates
                       tmp_dv_out.csv, tmp_dv_outf.
                4. APPENDING PREVIOUS DAY DECI-VIEW: appends to outfile.
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['day_loop_for_vis',]
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

def day_loop_for_vis(jdbeg,jdend,outfile,xref,comb,nifiles,infile_h,infile_t,class_lst,*,lyyyyjjj=True,tzone=None):
    # Include functions
    from CAMxtools.combine.combine import combine
    from CAMxtools.vis.calc_dv import calc_dv
    from CAMxtools.avg.get_davg_at_cells import get_davg_at_cells
    import numpy as np

    # Begin the script
    jdate = jdbeg
    if Path(outfile).exists(): os.remove(outfile)
    fout = open(outfile,'a')
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
      tracernames, indata2 = combine(None,comb,nifiles,infile,lverbose=False,loutf=False)
      jdatem1 = int((datetime.datetime.strptime(str(jdate),"%Y%j") + datetime.timedelta(days=-1)).strftime("%Y%j"))
      tmp_davg_outf = "tmp_davg_out.csv"
      tmp_dv_outf = "tmp_dv_out.csv"
      if (jdate > jdbeg):
        print ("  2. PROCESSING DAILY AVERAGE FOR PREVIOUS DAY")
        indata = np.append(indata1,indata2,axis=1)
        get_davg_at_cells(tracernames,indata,tmp_davg_outf,xref,jdatem1,tzone=tzone)
        print ("  3. PROCESSING DECI-VIEW FOR PREVIOUS DAY")
        calc_dv(tmp_davg_outf,tmp_dv_outf,class_lst)
        os.remove(tmp_davg_outf)
        print ("  4. APPENDING PREVIOUS DAY DECI-VIEW")
        fin = open(tmp_dv_outf,'r')
        header = fin.readline() # no. of header lines is one.
        if jdate == (jdbeg+1):
          fout.write('%s' % (str(header)))
        for line in fin:
          fout.write('%s' % (str(line)))
        fin.close()
        os.remove(tmp_dv_outf)
      indata1 = indata2
      jdate = int((datetime.datetime.strptime(str(jdate),"%Y%j") + datetime.timedelta(days=1)).strftime("%Y%j"))
    del indata2
    del indata1
    del indata
    fout.close()
