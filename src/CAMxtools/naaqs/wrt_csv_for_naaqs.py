__doc__ = r"""
.. _dumper
:mod:`dumper` -- CAMxtools dump module
============================================

... module:: wrt_csv_for_naaqs
    :platform: Unix, Windows
    :synopsis: Output to a csv file for naaqs
    :history:  Original version - zliu
               Heavily modified for CAMxtools - jjung 2/6/2018
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>
"""

__all__=['wrt_csv_for_naaqs',]

def wrt_csv_for_naaqs ( out_file,tracernames,data2csv,jdays,hours ) :
    # Include modules

    # set nx, ny, and nd
    nx = data2csv.shape[2]
    ny = data2csv.shape[1]
    ntracers = data2csv.shape[0]

    # write the ASCII file header
    fout = open( out_file , "w")
    fout.write('ICELL,JCELL,YJJJ,HR,IJCELL,')
    fout.write('%s' % (','.join(tracernames[nn] for nn in range (len(tracernames)) )))
    fout.write('\n')

    # write the ASCII file data
    for i in range ( nx ) :
      for j in range ( ny )  :
        ijcell = (i+1)*1000+j+1
        data1d = data2csv[:,j,i].flatten()
        fout.write('%d,%d,%d,%d,%d,' % (i+1,j+1,jdays[j,i],hours[j,i],ijcell) )
        fout.write('%s' % (','.join(str(data1d[ss]) for ss in range(ntracers))))
        fout.write('\n')
    fout.close()
    print('*** SUCCESS writing csv file')
    return
