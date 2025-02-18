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
domain= yantra.Domain3D((0,-dx,0),(l,(0.5*l)+dx,0.5*l),dx, grid_type = 'midway',periodicity={'x':1,'y':0,'z':1})
domain.nodetype[:,0,:]=1
domain.nodetype[:,-1,:]=1
#%% setup physics
g =0.0001 #acceleration like gravity
visc = 1./6.
poros =0.8
rho=1
K=1*poros**3/(150*(1-poros)**2)
domain_params={}
#domain_params['Fv']=[G,0]
domain_params['visc']=visc
domain_params['rho0']=rho
domain_params['poros']=poros
bc_params = {}

solver_params={'lattice':'D3Q19','collision_model':'srt','forcing_model':'mlga','magic_para':3./16.,'tauref':1}
ns = yantra.VolumeAveragedNavierStokes(domain,domain_params,bc_params,solver_params)
#%% run model
while ns.iters <=1000:
    Fv = np.zeros(ns.vector_shape)
    u = ns.u
    Fv[0] =  poros*rho*g - (rho*visc/K)*u[0]*poros
    Fv[1] =  poros*rho*0.- (rho*visc/K)*u[1]*poros    
    Fv[2] =  poros*rho*0 - (rho*visc/K)*u[2]*poros
    ns.Fv = Fv
    ns.advance()

#%%post process data
x,y,z = ns.meshgrid()

plt.figure()
plt.contourf(x[0,:,:],y[0,:,:],ns.u_mag[0,:,:])
plt.axis('image')
plt.title('velocity vector')
plt.colorbar()
plt.show()

plt.figure()
plt.quiver(x[0,:,:],y[0,:,:],ns.u[0][0,:,:],ns.u[1][0,:,:])
plt.axis('image')
plt.title('velocity')
plt.show()

x=x[0,0,:]
y=y[0,:,0]
u = analytic(g,rho,visc,poros,K,x,y,l,l/2)
plt.figure()
plt.title('velocity comparison')
plt.plot(y,ns.u[0][0,:,0],y,u,'o')
plt.legend(['LB','Analytical'])
plt.xlim([0,l/2])
plt.ylim([0,np.max(ns.u)+0.1*np.max(ns.u)])
plt.xlabel('distance')
plt.ylabel('velocity')
plt.show()
