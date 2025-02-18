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
Benchmark 19: steady state 2D heat transfer with conduction NAFEMS benchmark 
taken from comsol manual.The temprature at (0.04, 0.04) should be around 333 K
"""
#%% import library
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import LinearNDInterpolator
#%%
npts = 20 
dx= 0.08/npts
domain = yantra.Domain2D((-dx,0),(0.08,0.14),dx,grid_type='midway')
domain.nodetype[domain.y>=0.1,0]=1
domain.nodetype[domain.y<=0.04,0]=1
#%%
#domain params 
domain_params={}
domain_params['u']= (0,0)
domain_params['kappa'] =52 #W/(m.K)
domain_params['T']=273.15 #K
#set boundary conditions
#concentration gradient along x-axis all other boundary periodic
bc_params ={}
bc_params['left']=['tot_flux' ,5E5] #W/m2
bc_params['right']=['T',273.15]#K 
bc_params['top']=['T',273.15]#K
bc_params['bottom']=['T',273.15]#K
#solver parameters
solver_params={}
solver_params['interp']=0
ccd= yantra.ConvectionConduction(domain,domain_params,bc_params,solver_params)
#%%run 
ccd.steady_state(tol=1e-7)
#%%post process
x,y = ccd.meshgrid()
plt.figure()
plt.contourf(x,y,ccd.T,cmap='hot')
plt.colorbar()
plt.axis('image')
plt.figure()
plt.contour(x,y,ccd.T,cmap='hot')
plt.colorbar()
plt.axis('image')
#%%
x= x.flatten() 
y= y.flatten() 
T= ccd.T.flatten() 
InterpFunc=LinearNDInterpolator(np.array([list(x), list(y)]).T, list(T))
print('+++T @ (0.04,0.04):%s'%InterpFunc(0.04,0.04))
