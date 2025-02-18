#!/usr/bin/python
# -*- coding: utf-8 -*-
#=======================================================================================
#This File is part of Yantra: A lattice Boltzmann method based tool for multiscale/
#multiphyics simulations
#=======================================================================================
#
#Copyright (C) 2016-2017  <Author> Ravi A. Patel <Email> ravee.a.patel@gmail.com
#
#This program is free software: you can redistribute it and/or modify it under the
#terms of the GNU General Public License as published by the Free Software
#Foundation, either version 3 of the License, or any later version.
#This program is distributed in the hope that it will be useful, but WITHOUT ANY
#WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#problem defination and phreeqc results contributed by Hossein Fazeli (UiO)
#=======================================================================================

__doc__="""
Benchmark 12: Single species reaction kinetics using Yantra Phreeqc coupling.
"""
#%%import modules
import sys
import os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import matplotlib.pylab as plt
#%% read phreeqc results
def read_phreeqc():
    out=np.loadtxt('ext_results/benchmark12.csv',delimiter=',',skiprows=1)
    return out
#%%set domain
domain = yantra.Domain2D((0,0),(0.2,0.2/2),0.2/50.,grid_type='midway')
#%%set reactive transport model
#domain params
ux=5e-6
domain_params={}
domain_params['D']=1e-8 
domain_params['u']=[ux,0]
domain_params['database']='phreeqc.dat'
domain_params['phrqc_input_file']='benchmark12.phrq'
domain_params['solution_labels']=100002
#bc params
bc_params={}
bc_params['solution_labels']={'left':100001}
bc_params['top']=['flux',0]
bc_params['bottom']=['flux',0]
bc_params['left']=['c',0]
bc_params['right']=['open',0]
#solver parameters
solver_params={}
solver_params['collision_model']='srt'
solver_params['phrqc_flags']={}
solver_params['phrqc_flags']['smart_run']=False

rt= yantra.PhrqcReactiveTransport('AdvectionDiffusion', domain,domain_params,bc_params,solver_params)
#%%adjust bc for setting cauchy
for name in rt.fluid.components:
    comp=getattr(rt.fluid,name)
    comp.bc={'left':['flux',comp.bc['left'][1]*ux]}
#%%run model
A=[]
time=[]
while rt.time <=(100000):
    rt.advance()
    if (rt.iters%2==0):
        print("time: %s"%rt.time)
        time.append(rt.time)
        comp=getattr(rt.fluid,'[A]')
        A.append(comp.c[int(comp.ny/2),-1])
rt.kill_phrqc_workers()
#%%plot_results
plt.figure()
#plot lb results
plt.plot(time,np.array(A)*1000)
phrqc=read_phreeqc()
plt.plot(phrqc[:,0],phrqc[:,1],'o')
plt.legend(['yantra','phreeqc'])
plt.xlabel('time (s)')
plt.ylabel('[A] (mM)')
plt.show()
