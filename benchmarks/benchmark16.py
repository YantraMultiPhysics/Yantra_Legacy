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
Benchmark 16: Poiseulle's flow by gravity through \
forcing term and zero velocity through bounce-back in porous media
"""
#%% import library
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import matplotlib.pylab as plt
#%% define analytical solution
def analytic(G,rho,visc,poros,K,x,y,l,H):
    r = np.sqrt(poros/K) 
    u = (G*K/visc)*(1-(np.cosh(r*(y-H/2)))/np.cosh(r*H/2))
    return u

#%%   setup domain
l = 20
dx = 1
#adding additional dx to have out of domain solid nodes in top and bottom
domain= yantra.Domain2D((0,-dx),(l,(0.5*l)+dx),dx, grid_type = 'midway',periodicity={'x':1,'y':0})
domain.nodetype[0,:]=1
domain.nodetype[-1,:]=1
#%% setup physics
g =0.0001 #acceleration like gravity
visc = 1./6.
poros =.8
rho=1
K=1*poros**3/(150*(1-poros)**2)
domain_params={}
#domain_params['Fv']=[G,0]
domain_params['visc']=visc
domain_params['rho0']=rho
domain_params['u']=[0,0]
domain_params['poros']=1
bc_params = {}

solver_params={'lattice':'D2Q9','collision_model':'trt','forcing_model':'sc','magic_para':3./16.,'tauref':1}
ns = yantra.IncompressibleNavierStokes(domain,domain_params,bc_params,solver_params)
#%% run model
while ns.iters <=8000:
     Fv = np.zeros(ns.vector_shape)
     u = ns._u
     Fv[0] = 1* (poros*rho*g - (ns.rho*visc/K)*u[0]*poros)
     Fv[1] = 1* (poros*rho*0. - (ns.rho*visc/K)*u[1]*poros)
     ns._Fv = Fv
     ns.advance()

#%%post process data
x,y = ns.meshgrid()

plt.figure()
plt.contourf(x,y,ns.u[0])
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
u = analytic(g,rho,visc,poros,K,x,y,l,l/2)
plt.figure()
plt.title('velocity comparison')
plt.plot(y,ns.u[0][:,0],y[1:-1],u[1:-1],'o')
plt.legend(['LB','Analytical'])
plt.xlim([0,l/2])
plt.xlabel('distance')
plt.ylabel('velocity')
plt.show()
