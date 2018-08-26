__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: wrt_csv_for_vis
    :platform: Unix, Windows
    :synopsis: It extracts trc3d(ntracers,ny,nx) at the location
               in xref and write to csv output file.
    :details:  See the description under the function.
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['wrt_csv_for_vis',]
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
import pandas as pd
from pandasql import sqldf
import netCDF4 as ncdf4
import os

def wrt_csv_for_vis(tracernames, trc3d, xref, yyyyjjj, outfile):
    """
    Generate csv file which has daily average values of species used in
    visibility analysis at the grid cells specified in xref file
    Arguments:
       tracernames - species names to write to csv file
       trc3d - daily average value arrays of species used in visibility
               analysis
       xref - cross reference file that specify grid cells belong to class
              1 or 2 areas
       yyyyjjj - Julian date of processing day
       outfile - output csv file which has daily average values at grid 
                 cells specified
    """
    ntracers = trc3d.shape[0]
    header_out = "ICELL JCELL YJJJ IJCELL SRC VAL".split()
    sites = pd.read_csv(xref, sep=',')
    data_out=[]
    for ijcell in sites.IJCELL.unique() :
        i = int(ijcell/1000)-1
        j = ijcell%1000-1
        trc1d = []
        for s in range(ntracers) :
            trc1 = trc3d[s,j,i]
            trc1d.append(trc1)
            if (tracernames[s] == 'SeaSalt') :
              data_out.append(tuple([i+1,j+1,int(yyyyjjj),ijcell,'Sea Salt',trc1d[s]]))
            else :
              data_out.append(tuple([i+1,j+1,int(yyyyjjj),ijcell,tracernames[s],trc1d[s]]))
    data = pd.DataFrame.from_records(data_out, columns=header_out) #run queries

    qjoinsites = """select data.YJJJ, sites.GROUP1,
        sites.GNAME, data.SRC, data.VAL, data.ICELL, data.JCELL
        from data
        join sites using(IJCELL)"""

    data_w_sites = sqldf(qjoinsites, locals())

    data_w_sites.to_csv(outfile, index=False)

    return

def main():
    # Arguments
    tracernames = str(sys.argv[1]).split() # pass an argument closing quotation such as "NO NO2". It will be recognized as one argument.
    # trc3d is np array. It cannot be passed through a shell script argument.
    xref = str(sys.argv[2])
    yyyyjjj = str(sys.argv[3])
    outfile = str(sys.argv[4])

    # Extract trc3d(ntracers,ny,nx) at the location in xref and write to csv
    wrt_csv_for_vis(tracernames, trc3d, xref, yyyyjjj, outfile)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
