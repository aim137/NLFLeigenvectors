#! python
from aim.floquet.aux_functions import get_bands
from aim.floquet.nl_databases import NLdblist
from matplotlib import pyplot as plt

KSbands=get_bands()
lista=['job02_1' ,'job02_2','job02_3' ,'job02_4' ,'job02_5' ,'job02_6' ,'job02_7' ,'job02_8' ,'job02_9', 'job02_10', 'job02_11', 'job02_12']

for _kpt in range(KSbands.shape[1]):
  mydb=NLdblist(lista,_kpt+1,4) #_kpt runs from 0, _kpt+1 runs from 1 to BZ
  opt_evecs = mydb.find_QE(max_fl_mode=1)
  opt_evecs.output_for_fortran(4,_kpt+1,4,file='fortran_evecs_maxfl_1.in')
 
