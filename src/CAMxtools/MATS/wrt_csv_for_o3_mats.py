__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: wrt_csv_for_o3_mats
    :platform: Unix, Windows
    :synopsis: Output to an O3 MATS input
    :history:  Original version - zliu 10/24/2016
               Heavily modified for CAMxtools - jjung 1/5/2018
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>
"""

__all__=['wrt_csv_for_o3_mats',]

def wrt_csv_for_o3_mats ( out_file,ozone_mda8,jdays,lats,lons,cnvfac ) :
    # Include modules
    import datetime

    # set nx, ny, and nd
    nx = ozone_mda8.shape[2]
    ny = ozone_mda8.shape[1]
    nd = ozone_mda8.shape[0]

    # write the ASCII file
    fout = open( out_file , "w")
    fout.write('Day'+'\n')
    fout.write('_ID,_TYPE,   LAT ,  LONG  ,    DATE  , O3')
    fout.write('\n')

    for i in range ( nx ) :
      for j in range ( ny )  :
        gind = (i+1)*1000+j+1
        lat = lats[j,i]; lon = lons[j,i]
        for t in range ( nd ) :
          yyyymmdd = int(datetime.datetime.strptime(str(jdays[t]),"%Y%j").strftime("%Y%m%d"))
          data = ozone_mda8 [ t , j , i ] * cnvfac
          fout.write('%d,%s,%f,%f,%s,%f' % (gind,'""',lat,lon,yyyymmdd,data) )
          fout.write('\n')
    fout.close()
