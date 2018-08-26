__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: calc_dv
    :platform: Unix, Windows
    :synopsis: Takes daily averaged PM species concentrations and NO2 from 
               model and calculates delta deci-view, DDV by the model using
               tables shown in FLAG(2010) document.
    :details: 1. Model's PM species and NO2 daily average values used in
                 visibility calculation must be in input file, CONF.
                 CM (ug/m3) : Coarse PM = PM10 - PM2.5 = CPRM + CCRS
                 EC (ug/m3) : Elemental Carbon = PEC
                 _NH4_2SO4 (ug/m3) : Ammonium sulfate. Instead of adding NH4,
                                     full neutralization of sulfate is assumed.
                                     = 1.375*PSO4
                 NH4NO3 (ug/m3) : Ammonium nitrate. As ammonium sulfate, full
                                  neutralization of nitrate is assumed.
                                  = 1.29*PNO3
                 OM (ug/m3) : Organic Matters = POA
                 Sea Salt (ug/m3) : = NA + PCL
                 Soil (ug/m3) : = FPRM + FCRS
                 NO2 (ppm) : = NO2
              2. The input file, CONF has data for only one day as follows,
                 -------------------------------------------------------------
                 YJJJ,GROUP1,GNAME,SRC,VAL,ICELL,JCELL
                 2013271,127,CI_Badlands NP,CM,1.2708363533,170,19
                 2013271,127,CI_Badlands NP,EC,0.0338107198477,170,19
                 2013271,127,CI_Badlands NP,NH42SO4,0.811749696732,170,19
                 2013271,127,CI_Badlands NP,(NH4)2SO4,0.811749696732,170,19
                 2013271,127,CI_Badlands NP,NH4NO3,0.0424444079399,170,19
                 2013271,127,CI_Badlands NP,NO2,0.000344746134942,170,19
                 2013271,127,CI_Badlands NP,OM,0.272528797388,170,19
                 2013271,127,CI_Badlands NP,Sea Salt,0.0429391115904,170,19
                 2013271,127,CI_Badlands NP,Soil,0.696000516415,170,19
                                        .
                                        .
                                        .
                 2013271,1489,CI_Yellowstone NP,NO2,8.57501945575e-05,9,76
                 2013271,1489,CI_Yellowstone NP,OM,0.254105418921,9,76
                 2013271,1489,CI_Yellowstone NP,Sea Salt,0.041093274951,9,76
                 2013271,1489,CI_Yellowstone NP,Soil,0.533379435539,9,76
                 -------------------------------------------------------------
              3. Class I or II LiST file (Cross Reference file), CLST maps 
                 Class I or II areas that are  of interested in your project to
                 the sites listed in the tables in FLAG(2010). Many of them
                 are IMPROVE sites. For example, CLST file is as follows,
                 -------------------------------------------------------------
                 CI_Arches,Arches NP
                 CI_Bandelier,Bandelier NM
                 CI_Black_Canyon,Black Canyon of the Gunnison NP
                 CI_Bosque,Bosque del Apache Wilderness
                 CI_Canyonlands,Canyonlands NP
                 CI_Capitol_Reef,Capitol Reef NP                                                                     CI_Eagles_Nest,Eagles Nest Wilderness
                 CI_Flat_Tops,Flat Tops Wilderness
                                        .
                                        .
                                        .
                 -------------------------------------------------------------
              4. Output file, OUTF reports one maximum delta deci-view out of
                 multiple grid cells in each Class I or II area and corresponding
                 variables used such as growth factors (fL, fS, fSS), extiction
                 coefficients by natural background (bext_n) (EC_N, NO2_N, ...),
                 extiction coefficients by source (model) (bext_s) (EC_S, NO2_S,
                 ...), and total extinction coefficients (bext_t) (EC_T, NO2_T,
                 ...). OUTF looks as follows,
                 -------------------------------------------------------------
                 DATE,ClassI_II,IJCELL,fL(RH),fS(RH),fSS(RH),DDV,(NH4)2SO4_N,CM_N,EC_N,NH4NO3_N,NO2_N,OM_N,Rayleigh_N,Sea Salt_N,Soil_N,Total_N,(NH4)2SO4_T,CM_T,EC_T,NH4NO3_T,NO2_T,OM_T,Rayleigh_T,Sea Salt_T,Soil_T,Total_T,(NH4)2SO4_S,CM_S,EC_S,NH4NO3_S,NO2_S,OM_S,Rayleigh_S,Sea Salt_S,Soil_S,Total_S,
                 2013271,CI_Badlands NP,123114,1.600000,1.890000,2.050000,0.000027,5.01495840e-01,1.55400000e+00,2.00000000e-01,4.55412000e-01,0.00000000e+00,1.73940000e+00,9.00000000e+00,3.48500000e-02,4.90000000e-01,1.39751578e+01,5.01500513e-01,1.55401041e+00,2.00000931e-01,4.55413292e-01,1.45518581e-05,1.73940136e+00,9.00000000e+00,3.48500000e-02,4.90004181e-01,1.39751952e+01,4.67251338e-06,1.04092498e-05,9.31361868e-07,1.29155891e-06,1.45518581e-05,1.36196611e-06,0.00000000e+00,0.00000000e+00,4.18070567e-06,3.73992138e-05,
                                        .
                                        .
                                        .
                 2013271,CI_Yellowstone NP,126091,1.660000,1.980000,2.170000,0.000009,5.25320640e-01,1.37400000e+00,2.00000000e-01,4.77057000e-01,0.00000000e+00,1.73940000e+00,8.00000000e+00,1.10670000e-01,5.00000000e-01,1.29264476e+01,5.25322485e-01,1.37400525e+00,2.00000442e-01,4.77057015e-01,1.65042477e-06,1.73940058e+00,8.00000000e+00,1.10670000e-01,5.00002104e-01,1.29264595e+01,1.84539116e-06,5.25135947e-06,4.41686545e-07,1.53928764e-08,1.65042477e-06,5.84537490e-07,0.00000000e+00,0.00000000e+00,2.10390613e-06,1.18926984e-05,
                 -------------------------------------------------------------
    :warning: 1. Tables6-9 in FLAG(2010) are formatted as csv, which are used
                 as lookup tables. Make sure the location of files are correctly
                 referred.
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>

"""

__all__=['calc_dv',]
from scipy.io import netcdf
import numpy as np
import math
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
import csv
from collections import defaultdict
import collections
from math import log
import datetime
import pkg_resources
import os

def makehash():
    return collections.defaultdict(makehash)

def retro_dictify(frame):
    d = {}
    for row in frame.values:
        here = d
        for elem in row[:-2]:
            if elem not in here:
                here[elem] = {}
            here = here[elem]
        here[row[-2]] = row[-1]
    return d

def calbext_all (cnc_ps4, cnc_pn3, cnc_poa, cnc_pec, cnc_pfs, \
                 cnc_pcc, cnc_nox, cnc_ss, cnc_rl, \
                 frh_l, frh_s, frh_ss) :
    """
    Calculate bext for each speices using the equations shown in FLAG(2010)
    Arguments:
       cnc_* - concentrations of species, ug/m3 for PM, ppm for NOx,
               unit of bext for cnc_rl (No calculation for Rayleigh)
       frh_* - growth factros from tables in FLAG(2010)
    """
    # extinction efficiencies
    eeso4_l = 4.8
    eeso4_s = 2.2
    eeno3_l = 5.1
    eeno3_s = 2.4
    eepoa_l = 6.1
    eepoa_s = 2.8
    eepec = 10.0
    eepfs = 1.0
    eepcc = 0.6
    eenox = 333.0  # Assuming NO2 in ppm
    eesss = 1.7
    # split large and small particles
    if (cnc_ps4 < 20.0) :
       xx_ps4_l = cnc_ps4/20.0*cnc_ps4
    else :
       xx_ps4_l = cnc_ps4
    xx_ps4_s = cnc_ps4 - xx_ps4_l
    if (cnc_pn3 < 20.0) :
       xx_pn3_l = cnc_pn3/20.0*cnc_pn3
    else :
       xx_pn3_l = cnc_pn3
    xx_pn3_s = cnc_pn3 - xx_pn3_l
    if (cnc_poa < 20.0) :
       xx_poa_l = cnc_poa/20.0*cnc_poa
    else :
       xx_poa_l = cnc_poa
    xx_poa_s = cnc_poa - xx_poa_l
    # calculation of bext
    bext_ps4 = frh_l*(eeso4_l*xx_ps4_l)+frh_s*(eeso4_s*xx_ps4_s)
    bext_pn3 = frh_l*(eeno3_l*xx_pn3_l)+frh_s*(eeno3_s*xx_pn3_s)
    bext_poa = eepoa_l*xx_poa_l+eepoa_s*xx_poa_s
    bext_pec = eepec*cnc_pec
    bext_pfs = eepfs*cnc_pfs
    bext_pcc = eepcc*cnc_pcc
    bext_nox = eenox*cnc_nox
    bext_ss  = frh_ss*eesss*cnc_ss
    bext_rl  = cnc_rl
    bext = bext_ps4 + bext_pn3 + bext_poa + bext_pec + bext_pfs + \
           bext_pcc + bext_nox + bext_ss + bext_rl
    return bext, bext_ps4, bext_pn3, bext_poa, bext_pec, bext_pfs, \
           bext_pcc, bext_nox, bext_ss, bext_rl

def calc_dv(CONF,OUTF,CLST,*, lverbose=False):
    """
    Create OUTF which has delta deci-view, DDV, and bext's
     - Input and output files have text format.
    Arguments:
       CONF - daily average concentration from model at every grid cells
       OUTF - Max DDV and corresponding bext's at each class I or II area
       CLST - Class I or II LiST file mapping Class I or II to IMPROVE
              sites in FLAG(2010)
       lverbose - if True, print more messages to screen
    """
    # Constants
    fso4 = 1.375
    fno3 = 1.290
    months = "Jan Feb Mar Apr May Jun \
              Jul Aug Sep Oct Nov Dec".split()
    
    # Matching List file
    if lverbose: print ('Reading a matching list file that relates user defined Class I/II areas to the Class I/II areas from the FLAG2010 documents.')
    clsI = pd.Series.from_csv(CLST, sep=',', index_col=0, header=None)
    del CLST
    
    # Lookup tables
    if lverbose: print ('Reading FLAG2010 Table 6-9.')
    try: # set RPATH from the package
      RPATH = pkg_resources.resource_filename('CAMxtools', 'vis/data/')
    except:
      try: # set RPATH from RE Novato repository
        RPATH = '/models/camx/postproc/python/vis/Table/' #Reference path
      except:
        print('If you did not install CAMxtools or not using this from RE Novato office,')
        print('do either of this. Program stops because it cannot find a dir that has data.')
    TBL6 = RPATH + 'FLAG2010_Table6.csv'
    TBL7 = RPATH + 'FLAG2010_Table7.csv'
    TBL8 = RPATH + 'FLAG2010_Table8.csv'
    TBL9 = RPATH + 'FLAG2010_Table9.csv'
    natr = pd.read_csv(TBL6, sep=',', index_col=0, skiprows=1)
    FLRH = pd.read_csv(TBL7, sep=',', index_col=0, skiprows=1)
    FSRH = pd.read_csv(TBL8, sep=',', index_col=0, skiprows=1)
    FSSRH = pd.read_csv(TBL9, sep=',', index_col=0, skiprows=1)
    
    del TBL6
    del TBL7
    del TBL8
    del TBL9
    
    # Read Conc data
    if lverbose: print ('Reading concentration file')
    df = pd.read_csv(CONF,usecols=['YJJJ','GNAME','SRC','VAL','ICELL','JCELL'])
    df['IJCELL'] = df['ICELL']*1000 + df['JCELL']
    df = df[['YJJJ', 'GNAME', 'IJCELL', 'SRC', 'VAL']]
    model=retro_dictify(df)
    
    del CONF
    
    # Calculate natural bext
    if lverbose: print ('Calculating natural BEXT')
    bext_n = makehash()
    for index, row in natr.iterrows():
        cnc_ps4 = natr['(NH4)2SO4'][index]
        cnc_pn3 = natr['NH4NO3'][index]
        cnc_poa = natr['OM'][index]
        cnc_pec = natr['EC'][index]
        cnc_pfs = natr['Soil'][index]
        cnc_pcc = natr['CM'][index]
        cnc_nox = 0.
        cnc_ss  = natr['Sea Salt'][index]
        cnc_rl  = natr['Rayleigh'][index]
        for mon in months:
           bext_n[index]['Total'][mon], \
           bext_n[index]['(NH4)2SO4'][mon], \
           bext_n[index]['NH4NO3'][mon], \
           bext_n[index]['OM'][mon], \
           bext_n[index]['EC'][mon], \
           bext_n[index]['Soil'][mon], \
           bext_n[index]['CM'][mon], \
           bext_n[index]['NO2'][mon], \
           bext_n[index]['Sea Salt'][mon], \
           bext_n[index]['Rayleigh'][mon] = \
           calbext_all(cnc_ps4,cnc_pn3,cnc_poa,cnc_pec,cnc_pfs, \
                       cnc_pcc,cnc_nox,cnc_ss,cnc_rl, \
                       FLRH[mon][index],FSRH[mon][index], \
                       FSSRH[mon][index])
    if lverbose: print ('Finished natural BEXT')
    
    # Calculate source bext and total bext
    if lverbose: print ('Calculating total BEXT and source BEXT')
    bext_t = makehash(); bext_s = makehash()
    # Print model
    for date in model :
        yyyymmdd = datetime.datetime.strptime(str(date),"%Y%j").strftime("%Y%m%d")
        mm = int((int(yyyymmdd)%10000)/100)
        mon = months[mm-1]
        for clI_II in model[date] :
           if lverbose: print (" at {}".format(clI_II))
           classI = clsI[clI_II].strip()
           for ijcl in model[date][clI_II] :
              cnc_ps4 = model[date][clI_II][ijcl]['_NH4_2SO4']+natr['(NH4)2SO4'][classI]
              cnc_pn3 = model[date][clI_II][ijcl]['NH4NO3']+natr['NH4NO3'][classI]
              cnc_poa = model[date][clI_II][ijcl]['OM']+natr['OM'][classI]
              cnc_pec = model[date][clI_II][ijcl]['EC']+natr['EC'][classI]
              cnc_pfs = model[date][clI_II][ijcl]['Soil']+natr['Soil'][classI]
              cnc_pcc = model[date][clI_II][ijcl]['CM']+natr['CM'][classI]
              cnc_nox = model[date][clI_II][ijcl]['NO2']
              cnc_ss  = model[date][clI_II][ijcl]['Sea Salt']+natr['Sea Salt'][classI]
              cnc_rl  = natr['Rayleigh'][classI]
              bext_t[date][clI_II][ijcl]['Total'], \
              bext_t[date][clI_II][ijcl]['(NH4)2SO4'], \
              bext_t[date][clI_II][ijcl]['NH4NO3'], \
              bext_t[date][clI_II][ijcl]['OM'], \
              bext_t[date][clI_II][ijcl]['EC'], \
              bext_t[date][clI_II][ijcl]['Soil'], \
              bext_t[date][clI_II][ijcl]['CM'], \
              bext_t[date][clI_II][ijcl]['NO2'], \
              bext_t[date][clI_II][ijcl]['Sea Salt'], \
              bext_t[date][clI_II][ijcl]['Rayleigh'] = \
              calbext_all(cnc_ps4,cnc_pn3,cnc_poa,cnc_pec,cnc_pfs, \
                          cnc_pcc,cnc_nox,cnc_ss,cnc_rl, \
                          FLRH[mon][classI],FSRH[mon][classI], \
                          FSSRH[mon][classI])
              spcs = ['Total', '(NH4)2SO4', 'NH4NO3', 'OM', 'EC', 'Soil', 'CM', 'NO2', 'Sea Salt', 'Rayleigh']
              for spc in spcs :
                 bext_s[date][clI_II][ijcl][spc] = bext_t[date][clI_II][ijcl][spc] - bext_n[classI][spc][mon]
    if lverbose: print ('Finished total BEXT and source BEXT')
    
    # Calculate (1) bext for each day and grid cell for modeled species
    #           (2) delta deciview = 10*ln(bextmax+back_bext)/back_bext)
    
    if lverbose: print ('Find max across cells and date')
    mbext = makehash(); mdate = makehash(); mijcl = makehash(); ddv = makehash()
    for date in model :
        yyyymmdd = datetime.datetime.strptime(str(date),"%Y%j").strftime("%Y%m%d")
        mm = int((int(yyyymmdd)%10000)/100)
        mon = months[mm-1]
        for clI_II in model[date] :
           classI = clsI[clI_II].strip()
           mbext[clI_II] = 0.
           for ijcl in model[date][clI_II] :
              if (bext_t[date][clI_II][ijcl]['Total'] >= mbext[clI_II]) :
                mbext[clI_II] = bext_t[date][clI_II][ijcl]['Total']
                mdate[clI_II] = date
                mijcl[clI_II] = ijcl
                ddv[clI_II] = 10*log(
                    bext_t[date][clI_II][ijcl]['Total']/ \
                    bext_n[classI]['Total'][mon])
    
    if lverbose: print ('Finish finding max across cells and date')
    
    # Write output file
    
    if lverbose: print ('Write outputs')
    fout = open(OUTF,"w") 
    fout.write ('DATE,ClassI_II,IJCELL,fL(RH),fS(RH),fSS(RH),DDV,')
    date = list(bext_t.keys())[0]
    clI_II = list(bext_t[date].keys())[0]
    classI = clsI[clI_II].strip()
    ijcl   = list(bext_t[date][clI_II].keys())[0]
    for spc in sorted(bext_n[classI].keys()) :
       fout.write (spc + '_N,')
    for spc in sorted(bext_t[date][clI_II][ijcl].keys()) :
       fout.write (spc + '_T,')
    for spc in sorted(bext_s[date][clI_II][ijcl].keys()) :
       fout.write (spc + '_S,')
    fout.write ('\n')
    for clI_II in sorted(ddv.keys()) :
       classI = clsI[clI_II].strip()
       date = mdate[clI_II]
       ijcl = mijcl[clI_II]
       yyyymmdd = datetime.datetime.strptime(str(date),"%Y%j").strftime("%Y%m%d")
       mm = int((int(yyyymmdd)%10000)/100)
       mon = months[mm-1]
       if lverbose: print ('Max BEXT for {} at {} ( {} ) ,{}'.format(clI_II, date,mon,ijcl))
       fout.write('%d,%s,%d,%f,%f,%f,%f,' % (date,clI_II,ijcl,FLRH[mon][classI],\
                  FSRH[mon][classI],FSSRH[mon][classI],ddv[clI_II]))
       for spc in sorted(bext_n[classI].keys()) :
          fout.write('%.8e,' % (bext_n[classI][spc][mon]))
       for spc in sorted(bext_t[date][clI_II][ijcl].keys()) :
         fout.write('%.8e,' % (bext_t[date][clI_II][ijcl][spc]))
       for spc in sorted(bext_s[date][clI_II][ijcl].keys()) :
         fout.write('%.8e,' % (bext_s[date][clI_II][ijcl][spc]))
       fout.write ('\n')
    if lverbose: print ('Finish writing outputs')
    if lverbose: print ('################# Normal Completion #################')

def main():
    OUTF = str(sys.argv[1])
    CONF = str(sys.argv[2])
    CLST = str(sys.argv[3])
    calc_dv(CONF, OUTF, CLST, lverbose = False)

if __name__ == '__main__':
    # For internal use (no CAMxtools package installed), set the package path.
    try :
      package_path = os.environ['PACKAGE_PATH']
      sys.path.append(package_path)
    except :
      print ("PACKAGE_PATH environment variable is not set.")

    # Main
    main()
