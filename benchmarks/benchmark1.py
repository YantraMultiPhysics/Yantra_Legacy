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


__doc__= """ 
Benchmark 1: Comparision of 2D implementation of Yantra's Advection diffusion equation \
 with 1D analytical solution for infinite domain with constant concentration  \
 boundary conditions 
"""
#%% import modules 
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
from numpy import exp
from scipy.special import erfc 
import matplotlib.pylab as plt
from yantra import timer
#%% Functions defining analytical solution of the problem
#Ref:
def analytical_soln(cb,x,ts,D, u):
    c = (cb/2)*(erfc((x - u*ts)/(4*D*ts)**0.5)+
         erfc((x+ u * ts)/(4*D *ts)**0.5)*exp(u*x/D))
    return c
#%%create yantra domain instance
domain = yantra.Domain2D((0,0),(200,10),1, grid_type = 'nodal',periodicity={'x':0,'y':1})
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
clock = timer()
ade.run(time=5000)    
print('time taken for simulation: %s s'%clock.time_lapsed)
#%%plot results   
plt.figure()
plt.contourf(ade.x,ade.y,ade.c)
plt.colorbar()
plt.title('Concentration field')
plt.show()
#comparision with analytical model
plt.figure()
D=np.mean(D)
cany = analytical_soln(cb,ade.x,ade.time,D,ux)
plt.plot(ade.x,cany, label = 'analytical')
plt.plot(ade.x,ade.c[int(ade.ny/2),:],'--', label = 'Yantra')
plt.xlabel('distance (m)')
plt.ylabel(r'concentration (mol/m$^3$)')
plt.legend()
err = np.sqrt(np.sum((cany-ade.c[int(ade.ny/2),:])**2)/np.sum(cany**2))
plt.title("Time: %s S           collision model: %s         Error: %s"%(ade.time,ade.collision_model,err))    
plt.show()
