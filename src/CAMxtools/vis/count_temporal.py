_ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: combine
    :platform: Unix, Windows
    :synopsis: Takes annual daily DDV data by concatenating daily outputs from
               calc_dv.py as input. Returns no. of days larger than 1 DDV or 
               0.5 DDV.
    :details: 1. Input csv file format can be found in "calc_dv.py".
              2. Output csv file format is,
                 -------------------------------------------------------------
                 Area,threshold,daynum
                 CII_Bighorn Canyon NRA,0.5,1
                 CI_Northern Cheyenne Indian Reservation,0.5,2
                 -------------------------------------------------------------
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>
"""

__all__=['count_temporal',]
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
import collections
import operator

def makehash():
    return collections.defaultdict(makehash)

def count_temporal(INPF,OUTF):
    """
    Count no. of days larger than 1 DDV or 0.5 DDV
    Arguments:
       INPF - annual daily DDV csv file
       OUTF - csv file with no. of days larger than 1 or 0.5 for each class I or II. If no. of days = 0, there is no report.
    """
    f = open(INPF,'r')
    fout = open(OUTF,'w')
    ddv = makehash(); cnt = makehash()
    line = f.readline()
    # Print a header line
    fout.write('Area,threshold,daynum\n')
    header = line.strip().split(',')
    icol = {}
    for (l,val) in enumerate(header):
        if val == 'DATE' : icol['DATE'] = l
        if val == 'ClassI_II' : icol['ClassI_II'] = l
        if val == 'DDV' : icol['DDV'] = l
    for line in f:
        cols = line.strip().split(',')
        ddv[cols[icol['ClassI_II']]][cols[icol['DATE']]] = float(cols[icol['DDV']].strip())
    del INPF
    
    for clI_II in sorted(ddv.keys()):
        cnt[clI_II]['0.5'] = 0
        cnt[clI_II]['1.0'] = 0
        for yjjj in sorted(ddv[clI_II].keys()):
            if ddv[clI_II][yjjj] > 0.5:
                cnt[clI_II]['0.5'] += 1  
            if ddv[clI_II][yjjj] > 1.0:
                cnt[clI_II]['1.0'] += 1  

    for clI_II in sorted(cnt.keys()):
        for thrs in sorted(cnt[clI_II].keys()):
            fout.write('%s,%s,%s\n' % (str(clI_II),str(thrs),str(cnt[clI_II][thrs])))

def main():
    OUTF = str(sys.argv[1])
    INPF = str(sys.argv[2])
    count_temporal(INPF,OUTF)

if __name__ == '__main__':
    main()
