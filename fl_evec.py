from netCDF4 import Dataset
import numpy as np
import cmath as c

HBAR_eVfs = 6.58211899E-1
HA2EV = 27.2113834

class NLeVec:

  def get_tVecs(self,kpt,band):
    """kpt: type int from 1 to BZ
       band: type int from 1 to E%nbf
    """
    ds=Dataset(self.jobdirs[0]+'/ndb.RT_V_bands_K_section')
    size=ds.dimensions['RT_nbands'].size

    list_of_evecs = []
    for jobdir in self.jobdirs:
      evec=np.zeros(size,dtype=complex)
      ds=Dataset(jobdir+'/ndb.RT_V_bands_K_section')
      
      for i in range(size):
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

  def centr_mode(self,shifted_mode):
    return shifted_mode - (self.max_fl_mode)

  def shift_mode(self,mode):
    return mode - (self.max_fl_mode)
 
  def get_frequency(self):
    ds=Dataset(self.jobdirs[0]+'/ndb.Nonlinear')
    freq = float(ds['EXTERNAL_FIELD1'][5])*HA2EV
    period = 2*c.pi / (freq / HBAR_eVfs)
    return freq,period

  def build_exp_matrix(self):
    """builds matrix of exponentials with times
       and size given by self.times and self.tot_fl_modes
    """
    matrix = np.zeros((self.tot_fl_modes,self.tot_fl_modes),dtype=complex)
    for i,t in enumerate(self.times):  # rows - time
     for j in range(self.tot_fl_modes): # columns
       eta = self.shift_mode(j)
       matrix[i,j]=c.exp(-1j * eta * (self.freq/HBAR_eVfs) * t)
     
    inverse = np.linalg.inv(matrix)
    return inverse

  def __init__(self,jobdirs,kpt,band,qe_eV=None):
    """jobdirs: type list
       kpt: type int from 1 to BZ
       band: type int from 1 to E%nbf
    """
    self.jobdirs = jobdirs
    self.tvecs = self.get_tVecs(kpt,band)
    self.times = self.get_times()
    if qe_eV is not None:
      self.pvecs = self.get_pVecs(qe_eV)
    self.max_fl_mode  = int((len(jobdirs)-1)/2)
    self.tot_fl_modes = self.max_fl_mode * 2 + 1
    self.freq, self.period = self.get_frequency()
    self.exp_mat_m1 = self.build_exp_matrix()

  def printout_data(self,output_file='test1'):
  
    with open(output_file+'_evecs','w') as f:
     f.write('Time   c1.real             c1.imag              c2.real               c2.imag\n')
     for i,v in enumerate(self.tvecs):
       f.write(f'{self.times[i]}  {v[0].real}  {v[0].imag}  {v[1].real}  {v[1].imag}  \n')

    with open(output_file+'_pvecs','w') as f:
     f.write('Time   c1.real             c1.imag              c2.real               c2.imag\n')
     for i,e in enumerate(self.pvecs):
       f.write(f'{self.times[i]}  {e[0].real}  {e[0].imag}  {e[1].real}  {e[1].imag}  \n')

  def calculate_fVecs(qe_eV):
    """this will calculate the Floquet vectors for a given qe_eV
       and return """
    pass





