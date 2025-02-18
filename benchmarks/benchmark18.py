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
Benchmark 16: time varying porosity flow
"""
#%% import library
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import matplotlib.pylab as plt
#%% define analytical solution
def analytic(poros):
    return  1./poros

def poros(x):
    poros = 1. - 0.1-0.1*np.sin(2*np.pi*x/1.99)
    return poros
#%%   setup domain
l = 2
dx = 0.04
#adding additional dx to have out of domain solid nodes in top and bottom
domain= yantra.Domain2D((0,0),(1.99,(0.5*l)),dx, grid_type = 'nodal',periodicity={'x':1,'y':1})
#%% setup physics
x,y = domain.meshgrid()
visc = 1./6.
domain_params={}
#domain_params['Fv']=[G,0]
domain_params['visc']=visc
domain_params['rho0']=1
domain_params['poros']=poros(x)
bc_params = {}

solver_params={'lattice':'D2Q9','collision_model':'trt','forcing_model':'guo','magic_para':3./16.,'tauref':1}
ns = yantra.VolumeAveragedNavierStokes(domain,domain_params,bc_params,solver_params)
#%% run model
while ns.time <=20:
    Fv = np.zeros(ns.vector_shape)
    u = ns.u
    Fv[0] =  ns.poros*0.1*10 -ns.poros**2*(1-ns.poros)*10*(u[0])
    Fv[1] =  -ns.poros**2*(1-ns.poros)*10*u[1]
    ns.Fv = Fv
    ns.advance()
    if ns.iters %1000==0:
        print(ns.time)

#%%post process data

plt.figure()
plt.contourf(x,y,ns.u[0])
plt.axis('image')
plt.title('velocity vector')
plt.colorbar()
plt.show()

u = analytic(ns.poros)
plt.figure()
plt.title('velocity comparison')
plt.plot(x[1,:],ns.u[0][1,:],x[1,:],u[1,:])
plt.legend(['LB','Analytical'])
plt.xlabel('distance')
plt.ylabel('velocity')
plt.show()

plt.figure()
plt.title('porosity')
plt.plot(x[1,:],ns.poros[1,:],'-')
plt.legend(['LB','Analytical'])
plt.xlabel('distance')
plt.ylabel('porosity')
plt.show()

plt.figure()
plt.title('velocity comparison')
plt.plot(x[1,:],ns.u_darcy[0][1,:])
plt.legend(['LB','Analytical'])
plt.xlabel('distance')
plt.ylabel('velocity')
plt.show()