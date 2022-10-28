from aim.floquet.fl_evec import NLdblist
from matplotlib import pyplot as plt

lista=['job02_1' ,'job02_2','job02_3' ,'job02_4' ,'job02_5' ,'job02_6' ,'job02_7' ,'job02_8' ,'job02_9', 'job02_10']

mydb=NLdblist(lista,6,4)

#evecs=mydb.run_NL2FL(qe_eV=-1.1881694832306466)
evecs=mydb.run_NL2FL(max_fl_mode=2,qe_eV=-1.2881694832306466)

if (False):
  listof_qe = []
  listof_err = []
  for i in range(-5,6):
    qe_eV=(1+i/10000000000)*-1.188170527
    err=mydb.run_NL2FL(qe_eV=qe_eV)
    listof_qe.append(qe_eV)
    listof_err.append(err)
  
  plt.plot(listof_qe,listof_err, marker='.',color='gray',linestyle="-",ms=12)
  plt.show()

if (True):
  fig=plt.figure(figsize=(10,14))
  X,Y=evecs.calc_all_times_pVecs()
  plt.plot(X,Y[:,1].real)
  plt.plot(evecs.times,evecs.NL_in[:,1].real ,label='NL_in' ,marker='o',linestyle='',ms=12)
  plt.plot(evecs.times,evecs.NL_out[:,1].real,label='NL_out',marker='.',linestyle='',ms=12)
  plt.show()
