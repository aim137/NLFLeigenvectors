from aim.floquet.constants import HBAR_eVfs,HA2EV
from netCDF4 import Dataset
import numpy as np
import cmath as c

##################################################################################
### AUXILIARY FUNCTIONS ###
##################################################################################

def build_exp_matrix(listof_times=None,t0_fs=None,t_step=None,n_steps=None,max_fl_mode=None,freq=None,l_inv=True):
  """builds matrix of exponentials with times
     and size given by self.times and self.tot_fl_modes
  """
  if max_fl_mode is None:
    raise ValueError("Missing total number of FL modes on input to build_exp_matrix")
  if listof_times is None and (t0_fs is None or t_step is None or n_steps is None):
    raise ValueError("Missing times on input to build_exp_matrix")

  if listof_times is None:
    listof_times = []
    for step in range(-int(n_steps*0.06),n_steps):
      t = t0_fs + t_step * step
      listof_times.append(t)

  tot_fl_modes = 2 * max_fl_mode + 1
  matrix = np.zeros((len(listof_times),tot_fl_modes),dtype=complex)
  for i,t in enumerate(listof_times):  # rows - time
   for j in range(tot_fl_modes): # columns - eta
     eta = j - max_fl_mode
     matrix[i,j]=c.exp(-1j * eta * (freq/HBAR_eVfs) * t)

  if l_inv:
    inverse = np.linalg.inv(matrix[:tot_fl_modes,:])
    return matrix, inverse
  else:
    arrof_time = np.array(listof_times)
    return matrix, arrof_time

def get_bands(kpt=None,band=None):
  ds=Dataset('SAVE/ns.db1')
  KSevalues=ds['EIGENVALUES'][:]
  vb,vb_kpt,vbM=0,0,-111
  for _band in range(KSevalues.shape[2]):
    for _kpt in range(KSevalues.shape[1]):
      for _spin in range(KSevalues.shape[0]):
        if KSevalues[_spin,_kpt,_band] < 0. and KSevalues[_spin,_kpt,_band] > vbM:
           vb,vb__kpt,vbM = _band,_kpt,KSevalues[_spin,_kpt,_band]

  if kpt is None and band is None:
    return (KSevalues - vbM)*HA2EV
  else:
    return (KSevalues[0,kpt-1,band-1] - vbM)*HA2EV
