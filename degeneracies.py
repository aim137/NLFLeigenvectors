import os
import numpy as np
import cmath as c
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
from aim.floquet.constants import HBAR_eVfs,HA2EV
from aim.floquet.aux_functions import build_exp_matrix

##################################################################################
### CLASS FLdegeneracy ###
##################################################################################

class FLdegeneracy:

  def __init__(self,file_H,file_evecs,f=None,max_eta=1):#,T,f,max_eta,kpt,tag):
   #self.period = T
   self.max_fl_mode = max_eta
   self.tot_fl_modes = 2 * max_eta + 1
   #self.kpt = kpt
   #self.tag = tag
   #self.dir = 'results-'+str(self.tag)
   self.H = self.read_H_file(file_H)
   if f is None: 
    self.freq = self.H[0,0].real - self.H[1,1].real
   else:
    self.freq = f
   self.remove_time_derivative()
   self.list_QEs,self.list_vecs = self.read_evecs_file(file_evecs)
   self.average_H = self.build_average_H()

  def read_H_file(self,file_H):
    with open(file_H,'r') as f:
      _a=f.read()
      _b=_a.split('\n')
      _b.remove('')
      _rows = []
      for _c in _b: # _c is a row
        _d = _c.split(')') # _d is a row as list
        _d.pop(-1)
        _row = []
        for _e in _d:
          _re,_im = _e.lstrip().lstrip('(').split(',')
          _c = float(_re) + 1j * float(_im)
          _row.append(_c)
        _rows.append(_row)
      Mat = np.zeros((len(_rows),len(_row)),dtype=complex) 
      for i,_row in enumerate(_rows):
        Mat[i,:] = np.array(_row,dtype=complex)
    return Mat

  def remove_time_derivative(self):
    check = self.H.shape[0]/self.tot_fl_modes - self.H.shape[0]//self.tot_fl_modes
    if (check != 0.0): ValueError('Wrong number of fl modes on initialization of FLdegeneragy class')
    self.bands = int(self.H.shape[0]/self.tot_fl_modes)
    
    i_fks = 0
    for _band in range(self.bands):
      for _mode in range(self.tot_fl_modes):
        _shifted_mode = _mode - self.max_fl_mode
        self.H[i_fks,i_fks] += _shifted_mode * self.freq
        i_fks += 1
        
    
  def read_evecs_file(self,file_evecs):
    list_of_QE   = []
    list_of_vecs = []
    with open(file_evecs,'r') as f:
      f.readline()
      for line in f.readlines():
        _line_as_string = line.strip('\n').lstrip().rstrip()
        if _line_as_string[0] != '(': #it is a quasi-energy
          _qe = float(_line_as_string)
          list_of_QE.append(_qe)
        else: #it is a vector
          _line_as_list = _line_as_string.split(")")
          _line_as_list.pop(-1)
          _vec_as_list = []
          for _comp in _line_as_list:
            _re,_im = _comp.lstrip().lstrip('(').split(',')
            _c = float(_re) + 1j * float(_im)
            _vec_as_list.append(_c)
          _vec_as_np = np.array(_vec_as_list,dtype=complex)
          list_of_vecs.append(_vec_as_np)
        
    return list_of_QE,list_of_vecs
 
  def build_average_H(self):
    dim = len(self.list_vecs)
    M = np.zeros((dim,dim),dtype=complex)
    for i in range(dim):
      for j in range(dim):
        M[i,j] = np.vdot(self.list_vecs[i],np.matmul(self.H,self.list_vecs[j]))              
    return M
  
  def diagonalize_average_H(self):
    vals,vecs=np.linalg.eig(self.average_H)
    return vals,vecs

  def find_band(self,fks_space=False):
    vals,vecs = self.diagonalize_average_H()
    chosen = abs(vals).argmin()
    coefficients = vecs[:,chosen]
    orig_vecs = np.array(self.list_vecs).transpose()
    new_vec = np.matmul(orig_vecs,coefficients)
    if fks_space:  return new_vec
    bands_modes_mat = new_vec.reshape((self.tot_fl_modes,self.bands),order='F')
    return bands_modes_mat


  def get_degen_vecs(self):
    list_reshaped = []
    for _vec in self.list_vecs:
      vec_reshaped = _vec.reshape((self.tot_fl_modes,self.bands),order='F')
      list_reshaped.append(vec_reshaped)
    return list_reshaped





