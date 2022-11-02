import os
import numpy as np
import cmath as c
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
from aim.floquet.constants import HBAR_eVfs,HA2EV
from aim.floquet.aux_functions import build_exp_matrix

##################################################################################
### CLASS FLeigenvectors ###
##################################################################################

class FLeigenvectors:

  def __init__(self,T,f,max_eta,listof_times,tag):
    self.period = T
    self.freq = f
    self.max_fl_mode = max_eta
    self.tot_fl_modes = 2 * max_eta + 1
    self.times = np.array(listof_times)
    self.tag = tag
    self.dir = 'figs-'+str(self.tag)

  def store_results(self,NL_in,NL_out,FL_vecs,FL_qe,err):
    self.NL_in = NL_in
    self.NL_out = NL_out
    self.FL_vecs = FL_vecs
    self.FL_qe = FL_qe
    self.err = err
    
  def calc_all_times_pVecs(self,t_step=0.0025,n_steps=669):
    M,alltimes = build_exp_matrix(
                 t0_fs=self.times[0],
                 t_step=t_step,
                 n_steps=n_steps,
                 max_fl_mode=self.max_fl_mode,
                 freq=self.freq,
                 l_inv=False)

    NL_out_alltimes = np.matmul(M,self.FL_vecs)
    return alltimes,NL_out_alltimes

  def plot_realtime(self,band_to_plot=1,t_step=0.0025,n_steps=669):
    X,Y=self.calc_all_times_pVecs(t_step=t_step,n_steps=n_steps)

    fig,axes=plt.subplots(2)
    fig.set_size_inches(8.3,5.8)
    fig.suptitle(f'Time-dependent projection over Kohn-Sham state: {band_to_plot+1}')

    axes[0].plot(X,Y[:,band_to_plot].real,label='FL_calculated',color='tab:blue')
    axes[0].plot(self.times,self.NL_in[:,band_to_plot].real ,label='NL_in' ,marker='o',linestyle='',color='tab:olive',ms=12)
    axes[0].plot(self.times,self.NL_out[:,band_to_plot].real,label='NL_out',marker='.',linestyle='',color='tab:blue',ms=11)
    axes[1].plot(X,Y[:,band_to_plot].imag,label='FL_calculated',color='tab:blue')
    axes[1].plot(self.times,self.NL_in[:,band_to_plot].imag ,label='NL_in' ,marker='o',linestyle='',color='tab:olive',ms=12)
    axes[1].plot(self.times,self.NL_out[:,band_to_plot].imag,label='NL_out',marker='.',linestyle='',color='tab:blue',ms=11)
    axes[1].set_xlabel('Time (fs)')
    axes[0].set_ylabel(f'Re[ d_{band_to_plot+1} ]')
    axes[1].set_ylabel(f'Im[ d_{band_to_plot+1} ]')
    os.system(f'if [ ! -d {self.dir} ]; then mkdir {self.dir};fi')
    plt.savefig(f'{self.dir}/fig-real_time_projection_over_KS_state_{band_to_plot+1}.pdf')
    plt.close()

  def plot_floquet(self):
    _fl_vec = self.FL_vecs.reshape((1,self.FL_vecs.size),order='F')
    fig = plt.figure()
    fig.set_size_inches(8.3,3.9)
    ax = fig.gca()
    fig.suptitle('Projection over Floquet-Kohn-Sham space')

    cax = ax.matshow(abs(_fl_vec), cmap='YlOrBr', norm=pltcol.LogNorm(vmin=1.e-10,vmax=1.))
    #cax = ax.matshow(abs(_fl_vec), cmap='gist_stern_r', vmin=1.e-16,vmax=1.)
    fig.colorbar(cax,orientation='horizontal')

    lof_labels=[]
    for band in range(self.FL_vecs.shape[1]):
      for i in range(self.FL_vecs.shape[0]):
        tic='        '+str(i-self.max_fl_mode)
        lof_labels.append(tic)
    ax.set_xticks([x-0.5 for x in range(_fl_vec.shape[1])])
    ax.set_xticklabels(lof_labels)
    ax.set_yticks([])
    ax.set_yticklabels([])

    ax.text(4,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode,0]):.1e}', fontsize=15)
    ax.text(0,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode+1,0]):.1e}', fontsize=15)
    ax.text(10,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode-1,1]):.1e}', fontsize=15)
    ax.text(13.5,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode+1,1]):.1e}', fontsize=15)
    plt.grid(which='major',lw=1.,color='black')

    os.system(f'if [ ! -d {self.dir} ]; then mkdir {self.dir};fi')
    plt.savefig(f'{self.dir}/fig-FKS_projection.pdf')
    plt.close()

  def output(self):
    os.system(f'if [ ! -d {self.dir} ]; then mkdir {self.dir};fi')
    with open(f'{self.dir}/output-{self.tag}.dat','w') as f:
      for k,v in self.__dict__.items():
        f.write('--------------------------\n')
        f.write(f'Attribute: {k}\n')
        f.write('Values:\n')
        f.write(str(v)+'\n')
