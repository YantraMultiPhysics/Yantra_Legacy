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
import matplotlib.pylab as plt
#%% define analytical solution
def analytic(G,rho,neu,x,y,lx,ly):
    y=y-(ly/2)
    u = G*ly**2/(8*neu*rho)*(1-4*(y**2/ly**2))
    return u

#%% setup domain
l = 20
dx = 1
#adding additional dx to have out of domain solid nodes in top-bottom and front-back
domain= yantra.Domain3D((0,-dx,-dx),(l,0.5*l+dx,0.5*l+dx),dx, grid_type = 'midway',periodicity={'x':1,'y':0,'z':1})
domain.nodetype[:,0,:]=1
domain.nodetype[:,-1,:]=1
#%% setup physics
G = 0.01
visc = 1./6.
domain_params={}
domain_params['Fv']=[G,0,0]
domain_params['visc']=visc
domain_params['rho0']=1
bc_params = {}
solver_params={'lattice':'D3Q19','collision_model':'trt','forcing_model':'guo','magic_para':3./16.,'tauref':1}
ns = yantra.IncompressibleNavierStokes(domain,domain_params,bc_params,solver_params)
#%% run model
ns.steady_state(tol=1e-6,check_iters=1,err_cal_method='l2norm')
#%%post process data
x,y,z = ns.meshgrid()
loc_z= int(ns.nz/2)
loc_y= int(ns.ny/2)
loc_x= int(ns.nx/2)
plt.figure()
umag=(ns.u[0]**2+ns.u[1]**2+ns.u[2]**2)**0.5
plt.contourf(x[loc_z,:,:],y[loc_z,:,:],umag[loc_z,:,:])
plt.axis('image')
plt.title('velocity vector')
plt.colorbar()
plt.show()

plt.figure()
plt.quiver(x[loc_z,:,:],y[loc_z,:,:],ns.u[0][loc_z,:,:],ns.u[1][loc_z,:,:])
plt.axis('image')
plt.title('velocity')
plt.show()

x=domain.x
y=domain.y#y[loc_z,:,loc_x]
u = analytic(G,ns.rhof,visc,x,y,l,l/2)
plt.figure()
plt.title('velocity comparison')
plt.plot(y,ns.u[0][loc_z,:,loc_x],'o',y,u)
plt.legend(['LB','Analytical'])
plt.xlabel('distance')
plt.ylabel('velocity')
plt.xlim([0,l/2])
plt.show()
print('error(in %):', ns.l2norm(ns.u[0][loc_z,1:-1,loc_x],u[1:-1])*100) #IGNORE BOUNDARIES IN ERROR COMPUTATION
