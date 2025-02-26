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
Benchmark 4: Comparision of 3D Yantra's Multilevel diffusion implementation with \
1D analytical solution for infinite domain with constant concentration boudary 
conditions
"""
#%% import modules 
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
from scipy.special import erfc 
import matplotlib.pylab as plt
#%% Functions defining analytical solution of the problem
#Ref:
def analytical_soln(cb,x,ts,D):
    c =  cb * erfc(x/np.sqrt(4*D*ts))   
    return c
#%%create yantra domain instance
domain = yantra.Domain3D((0,0,0),(200,20,20),1, grid_type = 'midway',periodicity={'x':0,'y':1,'z':1})
#%% Define physics
cb = 1
D=0.1
#domain params 
domain_params={}
domain_params['D0'] = D
#set boundary conditions
#concentration gradient along x-axis all other boundary periodic
bc_params ={}
bc_params['left']=['c' ,cb]
bc_params['right']=['c',0.0]
#solver parameters
solver_params={}
solver_params['lattice']='D3Q7'
solver_params['collision_model']='TRT'  
solver_params['Deref']= domain_params['D0']/10
solver_params['magic_para']= 1./4.

#create physics instance
diff = yantra.MultilevelDiffusion(domain,domain_params,bc_params,solver_params)
#%%run model
diff.run(time=10000)    
#%%plot results   
plt.figure()
plt.contourf(diff.x,diff.y,diff.c[int(diff.nz/2),:,:])
plt.colorbar()
plt.title('Concentration field')
plt.show()
#comparision with analytical model
plt.figure()
D=np.mean(D)
cany = analytical_soln(cb,diff.x,diff.time,D)
plt.plot(diff.x,cany, label = 'analytical')
plt.plot(diff.x,diff.c[int(diff.nz/2),int(diff.ny/2),:],'--', label = 'Yantra')
plt.xlabel('distance (m)')
plt.ylabel(r'concentration (mol/m$^3$)')
plt.legend()
err = np.sqrt(np.sum((cany-diff.c[int(diff.nz/2),int(diff.ny/2),:])**2)/np.sum(cany**2))
plt.title("Time: %s S           collision model: %s         Error: %s"%(diff.time,diff.collision_model,err))    

 
    
