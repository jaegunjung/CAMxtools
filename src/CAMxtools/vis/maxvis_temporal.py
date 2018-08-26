_ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: combine
    :platform: Unix, Windows
    :synopsis: Takes annual daily DDV data by concatenating daily outputs from 
               calc_dv.py as input. Returns Max, 8th, W20, B20, and Min DDV days.
    :details: 1. Input csv file format can be found in "calc_dv.py".
              2. Output csv file format is,
                 -------------------------------------------------------------
                 Rank,DATE,ClassI_II,IJCELL,fL(RH),fS(RH),fSS(RH),DDV,(NH4)2SO4_N,CM_N,EC_N,NH4NO3_N,NO2_N,OM_N,Rayleigh_N,Sea Salt_N,Soil_N,Total_N,(NH4)2SO4_T,CM_T,EC_T,NH4NO3_T,NO2_T,OM_T,Rayleigh_T,Sea Salt_T,Soil_T,Total_T,(NH4)2SO4_S,CM_S,EC_S,NH4NO3_S,NO2_S,OM_S,Rayleigh_S,Sea Salt_S,Soil_S,Total_S,
                 Max,2012279,CII_Absaroka-Beartooth Wilderness,38075,1.910000,2.360000,2.590000,0.022588,6.25902720e-01,1.15200000e+00,2.00000000e-01,5.68438500e-01,0.00000000e+00,1.73940000e+00,9.00000000e+00,4.40300000e-02,4.10000000e-01,1.37397712e+01,6.26121159e-01,1.16083342e+00,2.02867614e-01,5.84264723e-01,6.28986584e-04,1.74035891e+00,9.00000000e+00,4.40300000e-02,4.11737234e-01,1.37708420e+01,2.18439441e-04,8.83341953e-03,2.86761380e-03,1.58262230e-02,6.28986584e-04,9.58913715e-04,0.00000000e+00,0.00000000e+00,1.73723372e-03,3.10708297e-02,
                 8th,2013124,CII_Absaroka-Beartooth Wilderness,39070,1.940000,2.410000,2.650000,0.008767,6.39127200e-01,1.15200000e+00,2.00000000e-01,5.80455000e-01,0.00000000e+00,1.73940000e+00,9.00000000e+00,4.50500000e-02,4.10000000e-01,1.37660322e+01,6.39464432e-01,1.15589978e+00,2.00962376e-01,5.85904235e-01,9.00195558e-05,1.73983188e+00,9.00000000e+00,4.50500000e-02,4.10903858e-01,1.37781066e+01,3.37231775e-04,3.89978243e-03,9.62376289e-04,5.44923550e-03,9.00195558e-05,4.31878210e-04,0.00000000e+00,0.00000000e+00,9.03858163e-04,1.20743819e-02,
                 W20,2013094,CII_Absaroka-Beartooth Wilderness,37077,1.950000,2.420000,2.660000,0.000317,6.41785920e-01,1.15200000e+00,2.00000000e-01,5.82868500e-01,0.00000000e+00,1.73940000e+00,9.00000000e+00,4.52200000e-02,4.10000000e-01,1.37712744e+01,6.41789860e-01,1.15217515e+00,2.00076778e-01,5.82963695e-01,1.83141639e-05,1.73942478e+00,9.00000000e+00,4.52200000e-02,4.10042982e-01,1.37717116e+01,3.94026076e-06,1.75148511e-04,7.67783149e-05,9.51954992e-05,1.83141639e-05,2.47828914e-05,0.00000000e+00,0.00000000e+00,4.29821957e-05,4.37141837e-04,
                 B20,2013103,CII_Absaroka-Beartooth Wilderness,39075,1.950000,2.420000,2.660000,0.000000,6.41785920e-01,1.15200000e+00,2.00000000e-01,5.82868500e-01,0.00000000e+00,1.73940000e+00,9.00000000e+00,4.52200000e-02,4.10000000e-01,1.37712744e+01,6.41785924e-01,1.15200000e+00,2.00000003e-01,5.82868529e-01,4.21147085e-09,1.73940000e+00,9.00000000e+00,4.52200000e-02,4.10000002e-01,1.37712745e+01,4.41405523e-09,4.60884175e-09,3.25280108e-09,2.90619704e-08,4.21147085e-09,1.16983734e-09,0.00000000e+00,0.00000000e+00,2.33241249e-09,4.90513887e-08,
                 MIN,2013083,CII_Absaroka-Beartooth Wilderness,36072,2.030000,2.530000,2.870000,0.000000,6.70928160e-01,1.15200000e+00,2.00000000e-01,6.09340500e-01,0.00000000e+00,1.73940000e+00,9.00000000e+00,4.87900000e-02,4.10000000e-01,1.38304587e+01,6.70928160e-01,1.15200000e+00,2.00000000e-01,6.09340531e-01,2.05066249e-09,1.73940000e+00,9.00000000e+00,4.87900000e-02,4.10000000e-01,1.38304587e+01,2.01156203e-10,9.55845403e-10,4.87061447e-10,3.06664236e-08,2.05066249e-09,1.62277303e-10,0.00000000e+00,0.00000000e+00,2.71260681e-10,3.47946845e-08,
                  Max,2012361,CII_BENTON LAKE NATIONAL WILDLIFE REFUGE,15147,2.410000,3.120000,3.520000,0.037530,8.27066880e-01,9.30000000e-01,2.00000000e-01,7.51201500e-01,0.00000000e+00,1.73940000e+00,9.00000000e+00,5.98400000e-02,3.50000000e-01,1.38575084e+01,8.27230488e-01,9.32919885e-01,2.02651674e-01,7.94774378e-01,6.69025436e-04,1.74011847e+00,9.00000000e+00,5.98400000e-02,3.51410097e-01,1.39096140e+01,1.63607800e-04,2.91988542e-03,2.65167357e-03,4.35728780e-02,6.69025436e-04,7.18468240e-04,0.00000000e+00,0.00000000e+00,1.41009723e-03,5.21056357e-02,
                                         .
                                         .
                                         .
                 -------------------------------------------------------------
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['maxvis_temporal',]
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

def maxvis_temporal(INPF,OUTF):
    """
    Find Max, 8th, W20, B20, and Min DDV days from annual daily DDV data
    Arguments:
       INPF - annual daily DDV csv file
       OUTF - csv file with Max, 8th, W20, B20, and Min DDV at each class I or II
    """
    f = open(INPF,'r')
    fout = open(OUTF,'w')
    ddv = makehash(); dll = makehash()
    line = f.readline()
    # Print a header line
    fout.write('%s%s' % ('Rank,',str(line)))
    header = line.strip().split(',')
    icol = {}
    for (l,val) in enumerate(header):
        if val == 'DATE' : icol['DATE'] = l
        if val == 'ClassI_II' : icol['ClassI_II'] = l
        if val == 'DDV' : icol['DDV'] = l
    for line in f:
        cols = line.strip().split(',')
        ddv[cols[icol['ClassI_II']]][cols[icol['DATE']]] = float(cols[icol['DDV']].strip())
        dll[cols[icol['ClassI_II']]][cols[icol['DATE']]] = line
    del INPF
    
    for clI_II in sorted(ddv.keys()):
        sbyval = sorted(ddv[clI_II].items(), key=operator.itemgetter(1), reverse = True)
        
        ndays_m1 = len(dll[cols[icol['ClassI_II']]]) - 1
        if (ndays_m1 == 363) or (ndays_m1 == 364):
           fout.write('%s%s' % ('Max,', str(dll[clI_II][sbyval[0][0]])))
           fout.write('%s%s' % ('8th,', str(dll[clI_II][sbyval[7][0]])))
           fout.write('%s%s' % ('W20,', str(dll[clI_II][sbyval[72][0]])))
           fout.write('%s%s' % ('B20,', str(dll[clI_II][sbyval[292][0]])))
           fout.write('%s%s' % ('MIN,', str(dll[clI_II][sbyval[ndays_m1][0]])))
        else:
           print ('The no. of days must be either 364 or 365')
           print ('Your no. of days is {}'.format(ndays_m1 + 1))
           exit()

def main():
    OUTF = str(sys.argv[1])
    INPF = str(sys.argv[2])
    maxvis_temporal(INPF,OUTF)
if __name__ == '__main__':
    main()
