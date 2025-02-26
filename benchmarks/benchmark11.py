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
Benchmark 11: C-S-H dissolution benchmark as described in chapter 4 Patel (2016), PhD thesis, UGent
"""

#%% import modules
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import pickle
from copy import deepcopy
import matplotlib.pylab as plt
#%% helper functions
def err_norm(a,b,tol):
    en =[]
    for k in list(a.keys()):
         en.append((np.sum((a[k]-b[k])**2)/np.sum((a[k])**2))**0.5) 
    return  np.nanmax(en)/tol 
#%%generate mesh
domain = yantra.Domain2D((0,0),(10e-6,1e-6),0.1e-6, grid_type = 'nodal')

#%%define physics
domain_params={}
domain_params['D0']=2.2e-10
domain_params['database']='cemdata07.dat'
domain_params['solution_labels']=100005
domain_params['phrqc_input_file']='benchmark11.phrq'
domain_params['ss_names']={'Tob_jen_ss':['CSHQ_TobH','CSHQ_TobD','CSHQ_JenH','CSHQ_JenD']}
domain_params['solid_phases']={
 'CSHQ_TobH': {'c':0.1041,'mvol':55.30e-3,'type':'diffusive'},
 'CSHQ_TobD': {'c':2.5050,'mvol':47.95e-3,'type':'diffusive'},
 'CSHQ_JenH': {'c':2.1555,'mvol':75.63e-3,'type':'diffusive'},
 'CSHQ_JenD': {'c':3.2623,'mvol':80.58e-3,'type':'diffusive'},}
domain_params['app_tort_params']=[1.,1./3.]  
#bc params
bc_params={}
bc_params['left']=['conc',0]
bc_params['right']=['flux',0]
bc_params['solution_labels']={'left':100001}
#generate physics
solver_params={}
solver_params['phrqc_flags']={}
solver_params['phrqc_flags']['smart_run']=True
solver_params['collision_model']='trt'
solver_params['const_cphi']=1
solver_params['cphi']=0.6*1./3.
solver_params['phrqc_smart_run_tol']=1e-8
rt=yantra.PhrqcReactiveTransport('MultilevelAdvectionDiffusion',domain,domain_params,bc_params,solver_params)
#%%run model
c_old = deepcopy( rt.fluid.get_attr('c'))
#pid params
KP=0.075;KI=0.175;KD=0.01;fac=1;tol=5e-2
taumin = 1
taumax = 1000
en1= en2 = en = 1.
rejections = 0.
tauref_old= 1
tauref=1
#final time and timesteps to save results
Ts=10
tsave = [0.1,0.5,1,2,3,4,5,6,7,8,9,10]
while rt.time < Ts:
    rt.advance()
    #adaptive timestepping using PID scheme
    if True:
        c_new = deepcopy( rt.fluid.get_attr('c'))
        en2 = en1;en1 =  en;en = err_norm(c_new,c_old,tol)
        pidfact = fac*(en1/en)**KP * (1.0/en)**KI * (en1**2/(en*en2))**KD
        tauref = 0.5 + (tauref_old-0.5)*(pidfact)
        if (en >1):
            tauref = min(tauref,tauref_old)
        if (tauref > taumin):
            rejections += 1
        tauref = max(tauref,taumin)
        tauref = min(tauref,taumax)
        rt.fluid.set_attr('tauref',tauref,component_dict=False)
        tauref_old = tauref
        c_old = deepcopy(c_new)
    if rt.iters%10==0:
        print ('Time: %s, Time step: %s, tau_ref: %s, err:%s'%(rt.time,rt.dt,tauref,en))
    #save stuff at some time steps
    if round(rt.time,4) in tsave:
        #save results
        results={}
        tsave.remove(round(rt.time,4))
        results['Ca'] = rt.fluid.Ca.c
        results['Cl'] = rt.fluid.Cl.c
        results['Si'] = rt.fluid.Si.c
        results['CSHQ_TobH']= rt.solid.CSHQ_TobH.c
        results['CSHQ_TobD']= rt.solid.CSHQ_TobD.c
        results['CSHQ_JenH']= rt.solid.CSHQ_JenH.c
        results['CSHQ_JenD']= rt.solid.CSHQ_JenD.c       
        results['pH']= rt.phrqc.selected_output()['pH']
        results['poros']= rt.solid.poros 
        results['x']=domain.x
        results['time']=round(rt.time,4)
        fname = 'results_%s.pkl'%round(rt.time,4)
        with open(fname,'wb') as fwrite:
            pickle.dump(results,fwrite)
        del results      
rt.kill_phrqc_workers()          
#%%plot results
loc_y = int(domain.ny/2)
plt.rc('text', usetex=False)
f, axarr = plt.subplots(2,2, sharex=False,sharey=False)
#Ca aqueous
axarr[0,0].set_title('Ca (aqueous)')
axarr[0,0].plot(np.array(domain.x)*1e6,rt.fluid.Ca.c[loc_y,:]*1e3)
axarr[0,0].set_xlabel(r'Distance (\mu m)')
axarr[0,0].set_ylabel('Concentration (mmol/lit)')
#Ca Si
axarr[0,1].set_title('Si (aqueous)')
axarr[0,1].plot(np.array(domain.x)*1e6,rt.fluid.Si.c[loc_y,:]*1e3)
axarr[0,1].set_xlabel(r'Distance (\mu m)')
axarr[0,1].set_ylabel('Concentration (mmol/lit)')
#Ca solid
Ca_s = (0.8333333*rt.solid.CSHQ_TobD.c[loc_y,:] + 0.6666667*rt.solid.CSHQ_TobH.c[loc_y,:] + 
        1.3333333*rt.solid.CSHQ_JenH.c[loc_y,:] + 1.5*rt.solid.CSHQ_JenD.c[loc_y,:])
axarr[1,0].set_title('Si (solid)')
axarr[1,0].plot(np.array(domain.x)*1e6,Ca_s)
axarr[1,0].set_xlabel(r'Distance (\mu m)')
axarr[1,0].set_ylabel('Concentration (mol/lit)')
#Si solid
Si_s = (0.6666667*rt.solid.CSHQ_TobD.c[loc_y,:] + 1.0*rt.solid.CSHQ_TobH.c[loc_y,:] + 
        1.0*rt.solid.CSHQ_JenH.c[loc_y,:] + 0.6666667*rt.solid.CSHQ_JenD.c[loc_y,:])
axarr[1,1].set_title('Si (solid)')
axarr[1,1].plot(np.array(domain.x)*1e6,Si_s)
axarr[1,1].set_xlabel(r'Distance (\mu m)')
axarr[1,1].set_ylabel('Concentration (mol/lit)')
plt.show()        
        
