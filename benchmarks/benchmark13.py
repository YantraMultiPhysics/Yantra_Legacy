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
import time
import matplotlib.pylab as plt
#%% define analytical solution
def analytic(G,rho,neu,x,y,lx,ly):
    y=y-(ly/2)
    u = G*ly**2/(8*neu*rho)*(1-4*(y**2/ly**2))
    return u

#%%   setup domain
l = 12
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
t0=time.time()
ns.steady_state()
t1=time.time()
print ("time taken:",(t1-t0))
mlups= (domain.nx*domain.ny*ns.iters*1e-6)/(t1-t0)
print ("MLUPS: ",mlups)
#%%post process data
x,y = ns.meshgrid()

plt.figure()
plt.contourf(x,y,(ns.u[0]**2+ns.u[1]**2)**0.5)
plt.axis('image')
plt.title('velocity vector')
plt.colorbar()
plt.show()

plt.figure()
plt.quiver(x,y,ns.u[0],ns.u[1])
plt.axis('image')
plt.title('velocity')
plt.show()

x=x[0,:]
y=y[:,0]
u = analytic(G,ns.rhof,visc,x,y,l,l/2)
plt.figure()
plt.title('velocity comparison')
plt.plot(y,ns.u[0][:,0],'o',y,u)
plt.legend(['LB','Analytical'])
plt.xlim([0,l/2])
plt.xlabel('distance')
plt.ylabel('velocity')
plt.show()
print('error(in %):', ns.l2norm(ns.u[0][1:-1,0],u[1:-1])*100) #IGNORE BOUNDARIES IN ERROR COMPUTATION
