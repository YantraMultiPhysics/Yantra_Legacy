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
Benchmark 20: Transient 2D heat transfer with conduction  benchmark with FEMLAB
taken from FEMLAB manual.The temprature at (0.1, 0.3) should be around 186.5 c for 712 element
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
npts = 23 #chosen to have same elements
dx= 0.3/npts
domain = yantra.Domain2D((0,0),(0.3,0.4),dx,grid_type='midway')
#%%
#domain params 
domain_params={}
domain_params['u']= (0,0)
domain_params['kappa'] =52 #W/(m.c)
domain_params['T']=0 #K
domain_params['rho'] =7850 #kg/m3
domain_params['cp'] =460 #W/(m.c)

#set boundary conditions
#concentration gradient along x-axis all other boundary periodic
bc_params ={}
bc_params['left']=['tot_flux' ,0] #W/m2
bc_params['bottom']=['T',1000]#c 
bc_params['top']=['T',1000]# W/(m2·°C)
bc_params['right']=['T',1000]# W/(m2·°C)
#solver parameters
solver_params={}
ccd= yantra.ConvectionConduction(domain,domain_params,bc_params,solver_params)
#%%run 
ccd.run(time=190)
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
print('+++T @ (0.04,0.04):%s'%InterpFunc(0.1,0.3))

