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
Benchmark 22: Lid driven cavity thermal flow. Benchmark example taken from
S. Chen et al., A simple lattice Boltzmann model for conjugate heat transfer research, 
Int. J. Heat Mass Transfer (2016)
"""
#%% import library
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
import numpy as np
#%%NS problem
#%%domain
L = 1
domain_ns = yantra.Domain2D((0,-0.1),(L,L),dx=1/100.,grid_type='midway')
x,y = domain_ns.meshgrid()
domain_ns.nodetype [y<0]=1
#%% set navier-stokes physics
dparams_ns= {'visc':0.71}
bcparams_ns={'left':['u',[0,0]],
             'right':['u',[0,0]],
             'bottom':['u',[0,0]],
             'top':['u',[1,0]]}
sparams_ns = {'tauref':1}
ns= yantra.IncompressibleNavierStokes(domain_ns,dparams_ns,bcparams_ns,sparams_ns)
#%%thermal problem
#%%domain
domain_th = yantra.Domain2D((0,-0.1),(1,1),dx=1/100.,grid_type='midway')
x,y = domain_th.meshgrid()
domain_th.nodetype = -2. *  (y<0)
#%%
dparams_th= {}
lbda_ref =  1
lbda_sol =  1
cp_ref= 1#lbda_ref/1#0.71 * lbda_ref/visc
cp_sol=1#lbda_sol/1
dparams_th= {}
dparams_th['u']=ns.u
dparams_th['kappa']=np.where(domain_th.nodetype==-2,lbda_sol,lbda_ref)
bcparams_th = {'left':['grad_T',0],
               'right':['grad_T',0],
               'bottom':['T',1],
               'top':['T',0]}
sparams_th ={'collision_model':'TRT','tauref':1}
th = yantra.ConvectionConduction(domain_th,dparams_th,bcparams_th,sparams_th)
dtCoupl=max(ns.dt,th.dt)
print(dtCoupl,ns.dt,th.dt)
#%%run model
itr_coupl=0
t_coupl=0
while itr_coupl < 150000:
#coupling terms
    Fv = np.zeros(ns.vector_shape)
    Fv[1,:,:]=    th.T * (ns.nodetype<=0)
    ns.Fv = Fv
    ns.run(time=ns.time+dtCoupl)
    itr_coupl+=1
    t_coupl+=dtCoupl    
#coupling terms
    th.u = ns.u
    th.run(time=th.time+dtCoupl)
    if (itr_coupl%100==0):    print('Time:',t_coupl,'Iters:',itr_coupl)
#%%plot
u= th.u.copy()
x,y = domain_th.meshgrid()
u[:,y<0]=np.NaN
plt.streamplot(x,y,u[0],u[1],density=0.7)
plt.axis('image')
plt.show()

#%%plot
x,y = domain_th.meshgrid()
plt.contour(x,y,th.T,15,cmap='jet')
plt.colorbar()
plt.axis('image')
plt.show()
#%% nusselt number
dTdY = th.grad_T[1].flatten()
x,y = domain_th.meshgrid()
x = list(x.flatten())
y = list(y.flatten())
from scipy.interpolate import LinearNDInterpolator
InterpFunc=LinearNDInterpolator(np.array([x,y]).T, dTdY)
xy = np.array([list(th.x),[0]*len(th.x)]).T
gradw = -np.mean(InterpFunc(xy))
Nu=gradw
print(Nu)