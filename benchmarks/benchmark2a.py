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

__doc__="""
Benchmark 2: Comparision of 3D Yantra's multilevel advection-diffusion implementation with \
1D analytical solution for infinite domain with constant concentration boundary \
conditions
"""
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
from numpy import exp
from scipy.special import erfc 
import matplotlib.pylab as plt
#%% Functions defining analytical solution of the problem
#Ref:
def analytical_soln(cb,x,ts,D, u):
    c =  cb * erfc(x/np.sqrt(4*D*ts))   
    c = (cb/2)*(erfc((x - u*ts)/(4*D*ts)**0.5)+
         erfc((x+ u * ts)/(4*D *ts)**0.5)*exp(u*x/D))
    return c
#%%create yantra domain instance
domain = yantra.Domain3D((0,0,0),(200,20,20),1, grid_type = 'midway',periodicity={'x':0,'y':0,'z':1})
#%% Define physics
cb = 1
ux=0.01
D=0.1
#domain params 
domain_params={}
domain_params['u']= (ux,0,0)
domain_params['D0'] = D
#set boundary conditions
#concentration gradient along x-axis all other boundary periodic
bc_params ={}
bc_params['left']=['c' ,cb]
bc_params['right']=['c',0.0]
#solver parameters
solver_params={}
solver_params['lattice']='D3Q7'
solver_params['collision_model']='trt' 
solver_params['Deref']= domain_params['D0']/10
solver_params['magic_para']= 1/4.

#create physics instance
ade = yantra.MultilevelAdvectionDiffusion(domain,domain_params,bc_params,solver_params)
#%%run model
ade.run(time=5000)    
#%%plot results   
plt.figure()
plt.contourf(ade.x,ade.y,ade.c[int(ade.nz/2),:,:])
cbar=plt.colorbar()
cbar.set_label('concentration (mol/m$^3$)')
plt.title('Time: %s'%ade.time)
plt.show()
#comparision with analytical model
plt.figure()
D=np.mean(D)
cany = analytical_soln(cb,ade.x,ade.time,D,ux)
plt.plot(ade.x,cany, label = 'analytical')
plt.plot(ade.x,ade.c[int(ade.nz/2),int(ade.ny/2),:],'--', label = 'Yantra')
plt.xlabel('distance (m)')
plt.ylabel(r'concentration (mol/m$^3$)')
plt.legend()
err = np.sqrt(np.sum((cany-ade.c[int(ade.nz/2),int(ade.ny/2),:])**2)/np.sum(cany**2))
plt.title("Time: %s S           collision model: %s         Error: %s"%(ade.time,ade.collision_model,err))    

