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

import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
import numpy as np
__doc__="""
Benchmark 15: Permeability computation for array of BBC cylinder.

Ref: Jones, B. and Feng, Y. T., (2011), Permeability assessment of heterogeous porous media
using lattice Boltzmann method, Particles 2011, Barcelona, Spain.

Contributor: Jonas Bentz (university of Koblenz-landau)
"""
#%% parameters for simulation
l = 50
dx=0.5
pin =1e-3
pout = 0.
tauref = 1.0
visc = 1./6.
rho0 = 1
gradP = (pout-pin)/l
#%%we run simulation for set of porosity to plot the relationship
porosity_list = np.linspace(0.3,0.95,10)
r_list=np.sqrt((1-porosity_list)/(2*np.pi))*l
k_computed = []
for poros,r in zip(porosity_list,r_list):
    print("computing for r= %s and porosity= %s"%(r,poros))
    #setup domain
    domain= yantra.Domain2D((0,0),(l,l),dx, grid_type = 'midway',periodicity={'x':0,'y':1})
    domain.draw_circle((l/2,l/2),r,1)
    domain.draw_circle((0,0),r,1)
    domain.draw_circle((l,0),r,1)
    domain.draw_circle((0,l),r,1)
    domain.draw_circle((l,l),r,1)
    x,y = domain.meshgrid()
    #setup physics
    domain_params={}
    domain_params['visc']=visc
    domain_params['rho0']=rho0
    bc_params = {}
    bc_params['left']=['p',[pin,0]]
    bc_params['right']=['p',[pout,0]]
    bc_params['top']=['u',[0,0]]
    bc_params['bottom']=['u',[0,0]]
    solver_params={'lattice':'D2Q9'}
    solver_params['tauref']=tauref
    ns = yantra.IncompressibleNavierStokes(domain,domain_params,bc_params,solver_params)
    #run model
    ns.steady_state()
    #compute permeability
    u_x = np.average(ns.u[0])
    k = -u_x/(gradP)/(2*r)**2*(1*visc*ns.rho0)
    k_computed.append(k)
    #plot velocity vector and hydrostatic pressure
    plt.figure()
    u = (ns.u[0]**2+ns.u[1]**2)**0.5
    u[domain.nodetype==1]=np.NaN
    plt.contourf(x,y,u)
    plt.axis('image')
    plt.title('velocity vector')
    plt.colorbar()
    plt.show()
    plt.figure()
    p = ns.p.copy()
    p[domain.nodetype==1]=np.NaN
    plt.contourf(x,y,p)
    plt.axis('image')
    plt.title('hydrostatic pressure')
    plt.colorbar()
    plt.show()

#%%plot solution
plt.figure()
plt.semilogy(porosity_list,k_computed,'o', label="LB result")

porosity_list = np.linspace(0.2,0.95,20)
#analytical expression of Lee and Yang (1997)
k_analytical1 = ((porosity_list**3*(porosity_list-0.2146))/(31*(1-porosity_list)**1.3))
#analytical expression of Gebart (1991)
k_analytical2 = (4/(9*np.pi*6**0.5))*(np.sqrt((np.pi/(2*3**0.5))/(1-porosity_list))-1)**2.5

plt.semilogy(porosity_list,k_analytical2,'--', label='Gebart (1991)') 
plt.semilogy(porosity_list,k_analytical1, label='Lee and Yang (1997)') 
plt.xlabel('porosity')
plt.ylabel('$k/d^2$')
plt.legend()
plt.show()