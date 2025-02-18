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

#%% import library
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
#%%   setup domain
l = 500
domain= yantra.Domain3D((0,0,0),(l,0.5*l,0.5*l),1, grid_type = 'midway',periodicity={'x':0,'y':0,'z':1})
#%% setup physics
pin = 0.001
pout = 0.0
visc = 1./6.
domain_params={}
domain_params['visc']=visc
domain_params['rho0']=1
bc_params = {}
bc_params['left']=['p',[pin,0.0,0]]
bc_params['right']=['p',[pout,0.0,0]]
bc_params['top']=['u',[0,0,0]]
bc_params['bottom']=['u',[0,0,0]]
solver_params={'lattice':'D3Q19','collision_model':'trt','magic_para':3./16.}
#solver_params['interp']=0

ns = yantra.IncompressibleNavierStokes(domain,domain_params,bc_params,solver_params)
#%% run model
t= yantra.timer()
ns.run(iters=1000)
mlups = (ns.ncells*ns.iters)/t.time_lapsed*1e-6
print('MLUPS:%s'%mlups) 