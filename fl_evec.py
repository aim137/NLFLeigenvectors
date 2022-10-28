from netCDF4 import Dataset
import numpy as np
import cmath as c
from datetime import datetime

HBAR_eVfs = 6.58211899E-1
HA2EV = 27.2113834

##################################################################################
### AUX FUNCTIONS ###
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

##################################################################################
### CLASS NLdblist ###
##################################################################################

class NLdblist:

  def get_tVecs(self,kpt,band):
    """kpt: type int from 1 to BZ
       band: type int from 1 to E%nbf
    """
    ds=Dataset(self.jobdirs[0]+'/ndb.RT_V_bands_K_section')
    self.basis_size=int(ds.dimensions['RT_nbands'].size)

    list_of_evecs = []
    for jobdir in self.jobdirs:
      evec=np.zeros(self.basis_size,dtype=complex)
      ds=Dataset(jobdir+'/ndb.RT_V_bands_K_section')
      
      for i in range(self.basis_size):
       evec[i]=complex(ds['V_bands'][0,0,kpt-1,band-1,i,0],ds['V_bands'][0,0,kpt-1,band-1,i,1])
    
      list_of_evecs.append(evec)
    
    return list_of_evecs

  def get_times(self):
    """gets times in fs as list
    """
    list_of_times = []
    for jobdir in self.jobdirs:
      ds=Dataset(jobdir+'/ndb.RT_V_bands')
      time_fs=float(ds['RT_TIMEs_NOW'][1])*2.418884326505/100 #AUT2FS
      list_of_times.append(time_fs)
    
    return list_of_times
 
  def get_frequency(self):
    ds=Dataset(self.jobdirs[0]+'/ndb.Nonlinear')
    freq = float(ds['EXTERNAL_FIELD1'][5])*HA2EV
    period = 2*c.pi / (freq / HBAR_eVfs)
    return freq,period

  def __init__(self,jobdirs,kpt,band,qe_eV=None):
    """jobdirs: type list
       kpt: type int from 1 to BZ
       band: type int from 1 to E%nbf
    """
    self.jobdirs = jobdirs[:]
    self.tvecs = self.get_tVecs(kpt,band)
    self.times = self.get_times()
    self.freq, self.period = self.get_frequency()

  def printout_data(self,output_file='test1'):
  
    with open(output_file+'_evecs','w') as f:
     f.write('Time   c1.real             c1.imag              c2.real               c2.imag\n')
     for i,v in enumerate(self.tvecs):
       f.write(f'{self.times[i]}  {v[0].real}  {v[0].imag}  {v[1].real}  {v[1].imag}  \n')

    with open(output_file+'_pvecs','w') as f:
     f.write('Time   c1.real             c1.imag              c2.real               c2.imag\n')
     for i,e in enumerate(self.pvecs):
       f.write(f'{self.times[i]}  {e[0].real}  {e[0].imag}  {e[1].real}  {e[1].imag}  \n')

# Processing part

  def setup_fl_space(self,max_fl_mode=3):
    self.max_fl_mode  = max_fl_mode
    self.tot_fl_modes = self.max_fl_mode * 2 + 1
    if not (len(self.tvecs) > self.tot_fl_modes): raise ValueError("You need more time steps than total fl modes")
    self.exp_mat_long,self.exp_mat_m1 = self.get_exp_matrix()

  def get_exp_matrix(self):
    """uses variables of class to set up a call to the
       external build_exp_matrix
    """
    M, Mm1 = build_exp_matrix(
             listof_times=self.times,
             max_fl_mode=self.max_fl_mode,
             freq=self.freq)
    return M, Mm1


# def build_exp_matrix(self):
#   """builds matrix of exponentials with times
#      and size given by self.times and self.tot_fl_modes
#   """
#   matrix = np.zeros((len(self.times),self.tot_fl_modes),dtype=complex)
#   for i,t in enumerate(self.times):  # rows - time
#    for j in range(self.tot_fl_modes): # columns - eta
#      eta = self.shift_mode(j)
#      matrix[i,j]=c.exp(-1j * eta * (self.freq/HBAR_eVfs) * t)
#    
#   inverse = np.linalg.inv(matrix[:self.tot_fl_modes,:])
#   return matrix, inverse
  
  def centr_mode(self,shifted_mode):
    return shifted_mode - (self.max_fl_mode)

  def shift_mode(self,mode):
    return mode - (self.max_fl_mode)
 

  def calc_pVecs(self,qe_eV):
    """Returns the periodic part of the 
       Floquet basis functions
    """
    mat_of_pvecs = np.zeros((len(self.tvecs),self.basis_size),dtype=complex)
    for i,v in enumerate(self.tvecs):
      mat_of_pvecs[i,:] = c.exp(+1j * qe_eV * self.times[i] / HBAR_eVfs) * v

    return mat_of_pvecs


  def calc_fVecs(self,qe_eV,mat_of_pvecs=None):
    """this will calculate the Floquet vectors for a given qe_eV
       and return """
    if self.exp_mat_m1 is None: raise AttributeError("You need to initialize the FL space first")
    if mat_of_pvecs is None:
      mat_of_pvecs = self.calc_pVecs(qe_eV)

    mat_of_fvecs = np.matmul(self.exp_mat_m1,mat_of_pvecs[:self.tot_fl_modes,:])
 
    return mat_of_fvecs

  def recalc_pVecs_via_fl(self,qe_eV=None,mat_of_pvecs=None,mat_of_fvecs=None):
    """Function to calculate pVecs via the obtained fVecs at all times steps,
       including those not used to generate those fVecs, i.e., outside the
       first period considered. This allow to determine whether the qe used 
       truly makes the tVecs periodic
    """
    if qe_eV is None and mat_of_pvecs is None: raise ValueError("Provide either qe or pVecs")
    if qe_eV is None and mat_of_fvecs is None: raise ValueError("Provide either qe or fVecs")
    if mat_of_pvecs is None:
      mat_of_pvecs = self.calc_pVecs(qe_eV)
    if mat_of_fvecs is None:
      mat_of_fvecs = self.calc_fVecs(qe_eV,mat_of_pvecs=mat_of_pvecs)
    
    mat_of_recalc_pvecs = np.matmul(self.exp_mat_long,mat_of_fvecs)

    return mat_of_recalc_pvecs
 
    
  def run_NL2FL(self,qe_eV,max_fl_mode=3,tag=None,iter_num=None):
    """run
    """
    if qe_eV is None: raise ValueError("You need the qe as input")

    if tag is None:
      tag = datetime.today().strftime('%Y%m%d-%H.%M.%S')
      if iter_num is not None:
        tag += '_iteration_'+str(iter_num)
    
    fl_eigenvectors = FLeigenvectors(self.period,self.freq,max_fl_mode,self.times)

    self.setup_fl_space(max_fl_mode)
    NL_in = self.calc_pVecs(qe_eV=qe_eV)
    FL_out = self.calc_fVecs(qe_eV,mat_of_pvecs=NL_in)
    NL_out = self.recalc_pVecs_via_fl(mat_of_pvecs=NL_in,mat_of_fvecs=FL_out)
    err = np.sum(np.abs(NL_out - NL_in))
    
    fl_eigenvectors.store_results(NL_in,NL_out,FL_out,qe_eV,err)
    
    return fl_eigenvectors


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
