__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: wrt_csv_for_pm_mats
    :platform: Unix, Windows
    :synopsis: Output to an O3 MATS input
    :history:  Original version - zliu 10/24/2016
               Heavily modified for CAMxtools - jjung 1/5/2018
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>
"""

__all__=['wrt_csv_for_pm_mats',]

def wrt_csv_for_pm_mats ( out_file,davg_pm,jdays,lats,lons,tracernames,cnvfac ) :
    # Include modules
    import datetime

    # set nx, ny, and nd
    nx = davg_pm.shape[-1]
    ny = davg_pm.shape[-2]
    nspc = davg_pm.shape[-3]
    nd = davg_pm.shape[-4]

    # set the order of species
    idx_crstl = tracernames.index("CRUSTAL")
    idx_nh4   = tracernames.index("NH4")
    idx_so4   = tracernames.index("SO4")
    idx_ec    = tracernames.index("EC")
    idx_no3   = tracernames.index("NO3")
    idx_oc    = tracernames.index("OC")
    idx_pm25  = tracernames.index("PM25")
    idx_cm    = tracernames.index("CM")
    idx_salt  = tracernames.index("SALT")

    # write the ASCII file
    fout = open( out_file , "w")
    fout.write('DAY'+'\n')
    fout.write('_ID,_TYPE,LAT,LONG,DATE,CRUSTAL,NH4,SO4,EC,NO3,OC,PM25,CM,SALT')
    fout.write('\n')

    for i in range ( nx ) :
      for j in range ( ny )  :
        gind = (i+1)*1000+j+1
        lat = lats[j,i]; lon = lons[j,i]
        for t in range ( nd ) :
          yyyymmdd = int(datetime.datetime.strptime(str(jdays[t]),"%Y%j").strftime("%Y%m%d"))
          data = davg_pm [ t , : , j , i ] * cnvfac
          fout.write('%d,%s,%f,%f,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f' % (gind,'""',lat,lon,yyyymmdd,data[idx_crstl],data[idx_nh4],data[idx_so4],data[idx_ec],data[idx_no3],data[idx_oc],data[idx_pm25],data[idx_cm],data[idx_salt]) )
          fout.write('\n')
    fout.close()
