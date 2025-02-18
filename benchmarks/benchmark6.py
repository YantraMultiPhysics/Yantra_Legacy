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
Benchmark 6: validation of multilevel diffusion model for domain with varying diffusion \
coefficient and porosity  using Yantra's 2D implementation. The benchmark example refers to example of \
Perko and Patel (2014) Modern Phy. C Fig 2(right).
"""
#%%import modules
import sys
import os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt

#%% domain generation
domain=yantra.Domain2D((0,0),(0.1,0.1),0.1/100,grid_type='midway')
domain.draw_rect((0.03,0.02),(0.06,0.04),-2)
domain.draw_rect((0.08,0.02),(0.04,0.04),-3)
domain.draw_rect((0.03,0.07),(0.06,0.06),-4)

#%%set Multilevel advection diffusion model
Dm= 1e-12*(domain.nodetype==-2)+\
    5e-10*(domain.nodetype==-3)+\
    1e-10*(domain.nodetype==-4)+\
    1e-11*(domain.nodetype==-1)
#Dm = 1e-10
poros =.1*(domain.nodetype==-2)+\
       .1*(domain.nodetype==-3)+\
       .5*(domain.nodetype==-4)+\
       .5*(domain.nodetype==-1)
D0=1e-12
#domain parameters
domain_params={}
domain_params['D0']=D0
domain_params['app_tort']=Dm/D0
domain_params['poros']=poros
domain_params['c']=1
#boundary parameters
bc_params = {}
bc_params['top']=['c',0.0]
bc_params['left']=['flux',0.0]
bc_params['right']=['flux',0.0]
bc_params['bottom']=['flux',0.0]
#solver parameters
solver_params={}
solver_params['Deref']=D0
solver_params['tauref']=1
solver_params['cphi_fact']=1./3.
solver_params['collision_model']='trt'

#create instance
ade=yantra.MultilevelAdvectionDiffusion(domain,domain_params,bc_params,solver_params)
#%% run model
ts = 200*24*3600
print(ts/ade.dt)
ade.run(time=ts)   
#%%plot results
plt.figure() 
plt.contourf(ade.x,ade.y,ade.c,50,cmap='jet')
plt.axis('image')
plt.colorbar()
plt.title(r'Concentration (mol/m$^3$)')
plt.show()
plt.figure() 
CS=plt.contour(ade.x,ade.y,ade.c,50,cmap='jet')
plt.axis('image')
plt.title(r'Contour: Concentration (mol/m$^3$)')
plt.clabel(CS, fontsize=8)
plt.show()