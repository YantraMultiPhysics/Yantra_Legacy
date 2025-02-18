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
Benchmark 13: Poiseulle's flow by applying pressure gradient through \
forcing term and zero velocity through bounce-back
"""
#%% import library
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
#%%   setup domain
l = 2500
dx = 1
#adding additional dx to have out of domain solid nodes in top and bottom
domain= yantra.Domain2D((0,-dx),(l,(0.5*l)+dx),dx, grid_type = 'midway',periodicity={'x':1,'y':0})
domain.nodetype[0,:]=1
domain.nodetype[-1,:]=1
#%% setup physics
G = 0.0001
visc = 1./6.
domain_params={}
domain_params['Fv']=[G,0]
domain_params['visc']=visc
domain_params['rho0']=1
bc_params = {}
solver_params={'lattice':'D2Q9','collision_model':'srt','forcing_model':'guo','magic_para':3./16.,'tauref':0.7}
ns = yantra.IncompressibleNavierStokes(domain,domain_params,bc_params,solver_params)
#%% run model
t= yantra.timer()
ns.run(iters=1000)
mlups = (ns.ncells*ns.iters)/t.time_lapsed*1e-6
print('MLUPS:%s'%mlups) 