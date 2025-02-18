#!/usr/bin/python
# -*- coding: utf-8 -*-
#=======================================================================================
#This File is part of Yantra Multiphysics: A lattice Boltzmann method based tool for multiscale/
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
Benchmark 24: Ragyleigh bernard convection with conjugate heat transfer.
Benchmark is taken from Pan X., etal. J. of Comp. Phy. 369 (2018)"Efficient monolithic projection method for time-dependent conjugate heat transfer".
Benchmark description can be found in section 3.1. Result can be found in table 1. Current parameter setting is for Gr=1e7 and ks/kf=5 for which Nu = 9.01.
Value obtained by running this script is 8.823035172857601 i.e 2.05% error which can be reduced by changing grid size.
"""
#%% import library
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
import numpy as np
#%%domain size
L = 1
dwall = 0.2
dx  = 1/225
#%% Pr and Gr
Pr =0.7
Gr= 1e7
#%%Ra
Ra=Pr*Gr
#%%ks/kf
kskf= 5
#%%other params
cp_f= 1
cp_s= 1
rho_f = rho_s =1
kf = 1
ks = kskf*kf
alpha_f= kf/(cp_f*rho_f)
visc = Pr * alpha_f
beta = 1
Th = 1
Tc= 0
g= -Ra * (visc * alpha_f)/(beta*(Th-Tc)*L**3)
#%%NS problem
#%%domain
domain_ns = yantra.Domain2D((0,0),(L+dwall,L),dx=dx,grid_type='midway')
x,y = domain_ns.meshgrid()
domain_ns.nodetype [x>1]=1
#%% set navier-stokes physics
dparams_ns= {'visc':visc}
bcparams_ns={'left':['u',[0,0]],
             'right':['u',[0,0]],
             'bottom':['u',[0,0]],
             'top':['u',[0,0]]}
sparams_ns = {'tauref':0.56,'collision_model':'TRT'}
ns= yantra.IncompressibleNavierStokes(domain_ns,dparams_ns,bcparams_ns,sparams_ns)
#%%thermal problem
#%%domain
domain_th = yantra.Domain2D((0,0),(L+dwall,L),dx=dx,grid_type='midway')
x,y = domain_th.meshgrid()
domain_th.nodetype [x>1]=-2
#%%
dparams_th= {}
dparams_th['u']=ns.u
dparams_th['kappa']=np.where(domain_th.nodetype==-2,ks,kf)
bcparams_th = {'left':['T',Tc],
               'right':['T',Th],
               'bottom':['grad_T',0],
               'top':['grad_T',0]}
sparams_th ={'collision_model':'TRT','tauref':1}
th = yantra.ConvectionConduction(domain_th,dparams_th,bcparams_th,sparams_th)
#%% some plotting functions
def plotstreamlines(display=False):
    global ns,t_coupl,itr_coupl
    u= ns.u.copy()
    x,y = ns.meshgrid()
    u[:,x>1]=np.NaN
    plt.figure()
    plt.streamplot(x,y,u[0],u[1],density=1.5,color=ns.u_mag)
    plt.axis('image')
    plt.title("Time:%.5f s"%t_coupl)
    if display:
        plt.show()
    else:        
        plt.savefig('thRBConjugateHeat/streamline_itr_%s.svg'%itr_coupl,dpi=1200)
        plt.close()

def plotisotherms(display=False):
    global th,t_coupl,itr_coupl
    x,y = domain_th.meshgrid()
    plt.figure()
    plt.contour(x,y,th.T,15,cmap='jet')
    plt.axis('image')
    plt.colorbar()
    plt.title("Time:%.5f s"%t_coupl)
    if display:
        plt.show()
    else:        
        plt.savefig('thRBConjugateHeat/isotherm_itr_%s.svg'%itr_coupl,dpi=1200)
        plt.close()
#%%run model
itr_coupl=0
t_coupl=0
dtCoupl=max(ns.dt,th.dt)
print(dtCoupl,ns.dt,th.dt)
#for movies keep this on
if not (os.path.isdir('thRBConjugateHeat')):
    os.mkdir('thRBConjugateHeat')
while itr_coupl < 250000:
#coupling terms
    Fv = np.zeros(ns.vector_shape)
    Fv[1,:,:]=    -g*ns.rho*beta*(th.T-Tc)  #bossineuq's approximation
    ns.Fv = Fv
    ns.run(time=ns.time+dtCoupl)
    itr_coupl+=1
    t_coupl+=dtCoupl
#coupling terms
    th.u = ns.u
    th.run(time=th.time+dtCoupl)
    if (itr_coupl%100==0):    print('Time:',t_coupl,'Iters:',itr_coupl)
    #for movies keep this on
    if(itr_coupl%5000==0):
        plotstreamlines()
        plotisotherms()
#%%
plotstreamlines(display=True)
plotisotherms(display=True)       
#%% nusselt number
dTdX = th.tot_flux[0].flatten()
x,y = domain_th.meshgrid()
x = x.flatten()
y = y.flatten()
xy = np.array([list(x),list(y)]).T
from scipy.interpolate import griddata
xyOut = ([1]*len(th.y),th.y)
InterpLineData=griddata(xy, dTdX,xyOut,method='cubic')
Nu=gradw = -np.mean(InterpLineData)
print(Nu)