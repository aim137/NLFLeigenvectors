from netCDF4 import Dataset
import numpy as np
import cmath as c

HBAR_eVfs = 6.58211899E-1
HA2EV = 27.2113834

class NLdb:

  def get_tVecs(self,kpt,band):
    """kpt: type int from 1 to BZ
       band: type int from 1 to E%nbf
    """
    ds=Dataset(self.jobdirs[0]+'/ndb.RT_V_bands_K_section')
    self.basis_size=ds.dimensions['RT_nbands'].size

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
 
  def get_pVecs(self,qe_eV):
    """Returns the periodic part of the 
       Floquet basis functions
    """
    list_of_pvecs = []
    for i,v in enumerate(self.tvecs):
      period_vec = c.exp(+1j * qe_eV * self.times[i] / HBAR_eVfs) * v
      list_of_pvecs.append(period_vec)

    return list_of_pvecs

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
    self.jobdirs = jobdirs
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
    self.exp_mat_long,self.exp_mat_m1 = self.build_exp_matrix()

  def build_exp_matrix(self):
    """builds matrix of exponentials with times
       and size given by self.times and self.tot_fl_modes
    """
    matrix = np.zeros((len(self.times),self.tot_fl_modes),dtype=complex)
    for i,t in enumerate(self.times):  # rows - time
     for j in range(self.tot_fl_modes): # columns - eta
       eta = self.shift_mode(j)
       matrix[i,j]=c.exp(-1j * eta * (self.freq/HBAR_eVfs) * t)
     
    inverse = np.linalg.inv(matrix[:self.tot_fl_modes,:])
    return matrix, inverse
  
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
    absolute_error = np.sum(np.abs(mat_of_recalc_pvecs-mat_of_pvecs))

    return mat_of_pvecs, mat_of_recalc_pvecs, absolute_error

    
  def run_NL2FL(self,qe_eV,max_fl_mode=3):
    """run
    """
    if qe_eV is None: raise ValueError("You need the qe as input")
    
    self.setup_fl_space(max_fl_mode)
    A,B,err = self.recalc_pVecs_via_fl(qe_eV=qe_eV)
    
    return err


