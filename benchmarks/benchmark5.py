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
#=======================================================================================

__doc__="""
Benchmark 5: Comparison of different collision models for domain with varying diffusion \
coefficient using Yantra's 2D implementation. The benchmark example refers to example of \
Perko and Patel (2014) PRE paper Fig 5 and results in Fig 6.
"""
#%% import necessray
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
os.environ['OMP_NUM_THREADS']=str(2)
    
#%% generate domain instance
domain = yantra.Domain2D((0,0), (0.1,0.05),0.001, grid_type = 'midway',periodicity={'x':0,'y':1})
x,_ = domain.meshgrid()
Dlow  = 1e-11
Dhigh = 2000 * Dlow
D = Dlow * (x<=0.02) + Dhigh*(x>0.02)*(x<0.08) +  Dlow * (x>=0.08)
#%% define physics
#domain params
domain_params={}
domain_params['D']=D
domain_params['D0']=1e-11
domain_params['app_tort']=D/1e-11
domain_params['c']=0.1 
bc_params = {}
bc_params['left']=['c',0.]
bc_params['right']=['flux',0]       
solver_params={}         
solver_params['Dref']=Dlow/5
solver_params['Deref']=Dlow/5
solver_params['lattice']='D2Q5'
solver_params['collision_model']='srt'
ade_srt=yantra.AdvectionDiffusion(domain,domain_params,bc_params,solver_params)
solver_params['collision_model']='diff_vel'
ade_diff_vel=yantra.AdvectionDiffusion(domain,domain_params,bc_params,solver_params)
solver_params['collision_model']='trt'
ade_trt=yantra.AdvectionDiffusion(domain,domain_params,bc_params,solver_params)
#%% run models
tf=500*3600*24 #500 days
ade_srt.run(time=tf)
print('srt collision model finished. Number of iterations: %s'%ade_srt.iters)
ade_trt.run(time=tf)
print('trt collision model finished. Number of iterations: %s'%ade_trt.iters)
ade_diff_vel.run(time=tf)
print('diff_vel collision model finished. Number of iterations: %s'%ade_diff_vel.iters)
#%% plot data to generate fig 6 of the paper
y = int(ade_srt.ny/2)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 20
plt.plot(x[y,:],ade_srt.c[y,:],label=ade_srt.collision_model.upper())
plt.plot(x[y,:],ade_trt.c[y,:],'-.',label=ade_trt.collision_model.upper())
plt.plot(x[y,:],ade_diff_vel.c[y,:],'--',label=ade_diff_vel.collision_model.upper())
plt.legend(loc=4)
plt.grid()
plt.axis('image')
plt.xlabel('Distance (m)')
plt.ylabel(r'concentration (mol/m$^3$)')
plt.title('Time %s days'%(ade_srt.time/(3600*24)))
plt.show()