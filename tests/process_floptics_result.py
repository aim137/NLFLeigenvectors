import numpy as np
from aim.floquet.fl_eigenvecs import FLeigenvectors



def file2vector(file,max_fl_mode=4):
  lof_complex_values = []
  with open (file,'r') as f:
    f.readline()
    for line in f.readlines():
      line = line.strip('\n')
      if line != ' ':
        _lof_string_values = line.split(',')
        _lof_string_values.remove('')
        for value in _lof_string_values:
          _aux_list = value.split(' ')
          _re,_im = filter(None,_aux_list)
          _number = float(_re) + 1j* float(_im)
          lof_complex_values.append(_number)
  
  _vector = np.array(lof_complex_values)
  tot_fl_modes = max_fl_mode * 2 + 1
  bands = int(_vector.shape[0]/tot_fl_modes)
  _fl_vec = _vector.reshape((tot_fl_modes,bands),order='F')
  return _fl_vec


for i in range(19):
  _tag =f'FL_eigenvectors35.{i+1}it.dat'

  vec_array = file2vector(_tag,max_fl_mode=4)
  vecs = FLeigenvectors.FromArrayOnly(vec_array,f'kpt35.{i+1}it')
  vecs.plot_floquet(vecs,labels='+2')
  #print(abs(vecs.FL_vecs[6,1]))
