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
#
#=======================================================================================


#%% import modules 
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
from yantra import timer
#%%create yantra domain instance
l=2500
domain = yantra.Domain2D((0,0),(l,l),1, grid_type = 'nodal',periodicity={'x':0,'y':1})
#%% Define physics
cb = 1
ux=0.01
D=0.1
#domain params 
domain_params={}
domain_params['u']= (ux,0)
domain_params['D'] = D
#set boundary conditions
#concentration gradient along x-axis all other boundary periodic
bc_params ={}
bc_params['left']=['c' ,cb]
bc_params['right']=['c',0.0]
#solver parameters
solver_params={}
solver_params['lattice']='D2Q5'
solver_params['collision_model']='trt' #other options 'trt' and 'diff_vel'
solver_params['Dref']= domain_params['D']/10
solver_params['magic_para']= 1./6.

#create physics instance
ade = yantra.AdvectionDiffusion(domain,domain_params,bc_params,solver_params)
#%%run model
t = timer()
ade.run(time=1000)    
mlups = (ade.ncells*ade.iters)/t.time_lapsed*1e-6
print('MLUPS:%s'%mlups) 