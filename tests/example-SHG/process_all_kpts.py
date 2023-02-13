#! python
from aim.floquet.aux_functions import get_bands
from aim.floquet.nl_databases import NLdblist
from matplotlib import pyplot as plt

KSbands=get_bands()
lista=['job02_1' ,'job02_2','job02_3' ,'job02_4' ,'job02_5' ,'job02_6' ,'job02_7' ,'job02_8' ,'job02_9', 'job02_10', 'job02_11', 'job02_12']

with open('results_evec.dat','w') as f_evec:
 with open('results_eval.dat','w') as f_eval:
  for _kpt in range(KSbands.shape[1]):
    mydb=NLdblist(lista,_kpt+1,4) #_kpt runs from 0, _kpt+1 runs from 1 to BZ
    opt_evecs = mydb.find_QE(tag=f'results_at_kpt{_kpt+1}')
    opt_evecs.plot_realtime(band_to_plot=0,t_step=0.004,n_steps=1340)
    opt_evecs.plot_realtime(band_to_plot=1,t_step=0.004,n_steps=1340)
    opt_evecs.plot_floquet(labels='+2')
    opt_evecs.output()
    f_eval.write(f'Kpt {_kpt+1:2}: KS eval [eV] = {mydb.KSev:.6f} --- FL qe [eV] = {opt_evecs.FL_qe:.6f} -- It.: {opt_evecs.nr_it:2} -- err in period = {opt_evecs.err:.2e}\n')
    f_evec.write(f'Kpt {_kpt+1:2}: FL evec: (val,0) = {abs(opt_evecs.FL_vecs[4,0]):.3e} -- (cond,+2) = {abs(opt_evecs.FL_vecs[6,1]):.3e}\n')
  
  
