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

  def __init__(self,file_H,file_evecs,max_eta):#,T,f,max_eta,kpt,tag):
   #self.period = T
   self.max_fl_mode = max_eta
   self.tot_fl_modes = 2 * max_eta + 1
   #self.kpt = kpt
   #self.tag = tag
   #self.dir = 'results-'+str(self.tag)
   self.H = self.read_H_file(file_H)
   self.freq = self.H[0,0].real - self.H[1,1].real
   self.remove_time_derivative()
   self.evec1,self.evec2 = self.read_evecs_file(file_evecs)

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
          print(_c)
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
    with open(file_evecs,'r') as f:
      
    return evec1,evec2
