_ = r"""                                                                        
.. _dumper                                                                      
:mod:`dumper` -- CAMxtools dump module                                         
============================================                                    
                                                                                
... module:: _data2fin
    :platform: Unix, Windows                                                    
    :synopsis: Takes tracernames, data2sav, attr_in and returns
               fin class with attributes.
    :details:
    :warning:
    :history:
... moduleauthor:: Jaegun Jung <jjung@ramboll.com>                              
"""                                               

__all__=['fclass','_data2fin',]

class fclass:
    pass

def _data2fin(data2sav, tracernames, attr_in):
    # Declare fin class
    fin = fclass()
  
    # Copy attributes from attr_in dictionary
    for key in attr_in.keys():
      setattr(fin,key,attr_in[key])
    setattr(fin,'NVARS',len(tracernames))
    setattr(fin,'VAR-LIST',''.join(map(lambda x: '{:<16s}'.format(x),tracernames)))
  
    # Construct fin variables attributes
    fin.variables = {}
    i = 0
    for name in tracernames: # GRIDDED
      if attr_in['FTYPE'] == 1:
        fin.variables[name] = data2sav[i,:,:,:,:]
      else: # BOUNDARY CONDITION
        fin.variables[name] = data2sav[i,:,:,:]
      i+=1

    return fin
