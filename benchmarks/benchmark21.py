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
Benchmark 21: Transient 2D heat transfer with conduction  with varying cp and k.
Benchmark from  Karani and Huber (2015), PRE 91 023304, figure (9)
"""
#%% import library
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
from matplotlib import rc
#rc('text', usetex=True)
#%%
npts = 15 #chosen to have same elements
dx=1/npts
domain = yantra.Domain2D((0,0),(1,3),dx,grid_type='midway',periodicity={'x':1,'y':0})
domain.draw_rect((0.5,0.5),(1,1),-1)
domain.draw_rect((0.5,1.5),(1,1),-2)
domain.draw_rect((0.5,2.5),(1,1),-3)
x,y = domain.meshgrid()
plt.pcolor(x,y,domain.nodetype)
#%%
#domain params 
domain_params={}
domain_params['u']= (0.0,0.0)
kappa  = 1*(domain.nodetype==-1)+0.1*(domain.nodetype==-2)+1*(domain.nodetype==-3)
domain_params['kappa'] =kappa #W/(m.c)
domain_params['T']=0 #K
domain_params['rho'] =1 #kg/m3
cp =  1*(domain.nodetype==-1)+0.1/3.*(domain.nodetype==-2)+1*(domain.nodetype==-3)
domain_params['cp'] =cp #W/(m.c)

#set boundary conditions
#concentration gradient along x-axis all other boundary periodic
bc_params ={}
bc_params['top']=['T',0]#c 
bc_params['bottom']=['T',1]#c
#solver parameters
solver_params={'collision_model':'srt','tauref':1}
ccd= yantra.ConvectionConduction(domain,domain_params,bc_params,solver_params)
#%%run
plt.figure() 
for t in [0.1,0.5,1.0,2,10]:
    ccd.run(time=t)
    plt.plot(ccd.y,ccd.T[:,int(ccd.nx/2)],'-o',label= 't=%s'%t)
plt.legend()
plt.ylabel('$\Theta$')
plt.xlabel('y/H')
#%%post process
x,y = ccd.meshgrid()
plt.figure()
plt.contourf(x,y,ccd.T,cmap='hot')
plt.colorbar()

