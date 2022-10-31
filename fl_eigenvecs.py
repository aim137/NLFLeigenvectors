import numpy as np
import cmath as c
from datetime import datetime
from aim.floquet.constants import HBAR_eVfs,HA2EV
from aim.floquet.aux_functions import build_exp_matrix

##################################################################################
### CLASS FLeigenvectors ###
##################################################################################

class FLeigenvectors:

  def __init__(self,T,f,max_eta,listof_times):
    self.period = T
    self.freq = f
    self.max_fl_mode = max_eta
    self.tot_fl_modes = 2 * max_eta + 1
    self.times = np.array(listof_times)

  def store_results(self,NL_in,NL_out,FL_vecs,FL_qe,err):
    self.NL_in = NL_in
    self.NL_out = NL_out
    self.FL_vecs = FL_vecs
    self.FL_qe = FL_qe
    self.err = err
    
  def calc_all_times_pVecs(self,t_step=0.0025,n_steps=469):
    M,alltimes = build_exp_matrix(
                 t0_fs=self.times[0],
                 t_step=t_step,
                 n_steps=n_steps,
                 max_fl_mode=self.max_fl_mode,
                 freq=self.freq,
                 l_inv=False)

    NL_out_alltimes = np.matmul(M,self.FL_vecs)
    return alltimes,NL_out_alltimes

  def plot(self):
    pass

  def output(self):
    pass
