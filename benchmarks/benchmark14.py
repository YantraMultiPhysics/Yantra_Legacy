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
Benchmark 14: Poiseulle's flow by applying pressure gradient \
and zero velocity through boundary conditions
"""
#%% import library
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
#%% define analytical solution
def analytic(pin,pout,rho,neu,x,y,lx,ly):
    dp = pin-pout
    p = pin-dp*x/lx
    G = dp/lx
    y=y-(ly/2)
    u = G*ly**2/(8*neu*rho)*(1-4*(y**2/ly**2))
    return p,u
#%%   setup domain
l = 20
dx=1.0
domain= yantra.Domain2D((0,0),(l,0.5*l),dx, grid_type = 'nodal')
#%% setup physics
pin =0.001
pout = 0
tauref = 1.0
visc = 1./6.
domain_params={}
domain_params['visc']=visc
domain_params['rho0']=1
bc_params = {}
bc_params['left']=['p',[pin,0]]
bc_params['right']=['p',[pout,0]]
bc_params['top']=['u',[0,0]]
bc_params['bottom']=['u',[0,0]]
solver_params={'lattice':'D2Q9','collision_model':'trt','magic_para':3./16.}
solver_params['tauref']=tauref
ns = yantra.IncompressibleNavierStokes(domain,domain_params,bc_params,solver_params)
#%% run model
ns.steady_state(tol=1e-6)
#%%post process data
x,y = ns.meshgrid()

plt.figure()
plt.contourf(x,y,(ns.u[0]**2+ns.u[1]**2)**0.5)
plt.axis('image')
plt.title('velocity vector')
plt.colorbar()
plt.show()

plt.figure()
plt.contourf(x,y,ns.p)
plt.axis('image')
plt.title('hydrostatic pressure')
plt.colorbar()
plt.show()

plt.figure()
plt.quiver(x,y,ns.u[0],ns.u[1])
plt.axis('image')
plt.title('velocity')
plt.show()

x=x[0,:]
y=y[:,0]
loc_x,loc_y = int((ns.nx)/2),int((ns.ny)/2)
p,u = analytic(pin,pout,ns.rhof,visc,x,y,l,l/2)
plt.figure()
plt.title('velocity comparison')
plt.plot(y,ns.u[0][:,loc_x],'o',y,u)
plt.legend(['LB','Analytical'])
plt.xlabel('distance')
plt.ylabel('velocity')
plt.show()

plt.figure()
plt.title('pressure comparison')
plt.plot(x,ns.p[loc_y,:],'o',x,p)
plt.legend(['LB','Analytical'])
plt.xlabel('distance')
plt.ylabel('pressure')
plt.show()

print('error(in %):', ns.l2norm(ns.u[0][:,0],u)*100) 

