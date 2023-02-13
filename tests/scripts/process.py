#! python

from aim.floquet.nl_databases import NLdblist
from matplotlib import pyplot as plt

lista=['job02_1' ,'job02_2','job02_3' ,'job02_4' ,'job02_5' ,'job02_6' ,'job02_7' ,'job02_8' ,'job02_9', 'job02_10', 'job02_11', 'job02_12']

mydb=NLdblist(lista,6,4)


if (False):
  listof_qe = []
  listof_err = []
  for i in range(-5,6):
    qe_eV=(1+i/10000000000)*-1.188170527
    fl_evecs=mydb.run_NL2FL(qe_eV=qe_eV)
    listof_qe.append(qe_eV)
    listof_err.append(fl_evecs.err)
  
  plt.plot(listof_qe,listof_err, marker='.',color='gray',linestyle="-",ms=12)
  plt.show()

if (False):
  evecs=mydb.run_NL2FL(max_fl_mode=4,qe_eV=-1.1881694832306466)
  evecs.plot_realtime(band_to_plot=0)
  evecs.plot_realtime(band_to_plot=1)
  
if (False):
  evecs=mydb.run_NL2FL(max_fl_mode=4,qe_eV=0.)
  evecs.plot_realtime(band_to_plot=0)
  evecs.plot_realtime(band_to_plot=1)

if (False):
  evecs=mydb.run_NL2FL(max_fl_mode=4,qe_eV=-2.3)
  evecs.plot_realtime(band_to_plot=0)
  evecs.plot_realtime(band_to_plot=1)
 
if (False):
  optimized_evecs = mydb.find_QE(qe_eV=-2.3,qe_thrs=1e-7,err_thrs=1e-8,max_fl_mode=4)  
  print(f'QE= {optimized_evecs.FL_qe} eV with accuracy {optimized_evecs.nr_acc} eV after {optimized_evecs.nr_it} iterations - error in periodicity = {optimized_evecs.err}')
  optimized_evecs.plot_realtime(band_to_plot=0)
  optimized_evecs.plot_realtime(band_to_plot=1)

if (True):
  optimized_evecs = mydb.find_QE(qe_thrs=1e-7,err_thrs=1e-8,max_fl_mode=4)  
  print(f'QE= {optimized_evecs.FL_qe} eV with accuracy {optimized_evecs.nr_acc} eV after {optimized_evecs.nr_it} iterations - error in periodicity = {optimized_evecs.err}')
  optimized_evecs.plot_realtime(band_to_plot=0,t_step=0.004)
  optimized_evecs.plot_realtime(band_to_plot=1,t_step=0.004)
  optimized_evecs.plot_floquet()
  optimized_evecs.output()
  
