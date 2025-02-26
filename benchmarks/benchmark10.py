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
Benchmark 10: Plus shaped portlandite dissolution with geometry update. The results  \
of this benchmark are published in Patel et al (2014), Phy & Chem. of Earth.
"""

#%% import modules
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import time as timer
import matplotlib.pylab as plt

class MyRT(yantra.PhrqcReactiveTransport):
    def modifyeq_skip_phase(self):
        out = {}
        out['portlandite']= 0. * np.zeros(self.solid.shape)
        return out
    
if __name__=='__main__':
    #%%
    #generate mesh
    domain = yantra.Domain2D((0,0),(150e-6,150e-6),1e-6, grid_type = 'midway')
    #make plus geometry
    domain.draw_rect((75e-6,75e-6),(30e-6,7.5e-6),idx=1)
    domain.draw_rect((75e-6,75e-6),(7.5e-6,30e-6),idx=1)
    pqty = 5.*(domain.nodetype>0)
    #%%
    #domain params
    domain_params={}
    domain_params['D']=1e-9
    domain_params['database']='cemdata07.dat'
    domain_params['phrqc_input_file']='benchmark10.phrq'
    domain_params['solution_labels']=100001
    domain_params['eq_names']=['portlandite']
    domain_params['solid_phases']={'portlandite':{'type':'non_diffusive','mvol':1,'c':pqty}}
    domain_params['voxel_vol']=5
    #solver parameters
    solver_params={}
    solver_params['collision_model']='srt'
    solver_params['phrqc_flags']={}
    solver_params['phrqc_flags']['only_interface']=True
    solver_params['phrqc_flags']['smart_run']=True
    solver_params['phrqc_smart_run_tol']=1e-9
    #solver_params['phrqc_flags']['smart_run']=True
    #solver_params['phrqc_smart_run_tol']=1e-8
    rt= MyRT('AdvectionDiffusion', domain, domain_params,{},solver_params)
    #%%run model
    t0 = timer.time()
    time=[]
    AvgCa =[]
    while rt.time<=60:
        rt.advance()
        if (rt.iters%50==0):
            print ("+++Time: %s"%rt.time)
            print ("Active Phreeqc Nodes: %s"%rt.phrqc.nactive)
            print ("Solid nodes: %s"%rt.solid.nsolids)
            AvgCa.append(np.sum(rt.fluid.Ca.c)/np.sum(rt.fluid.Ca.nodetype<=0)) 
            time.append(rt.time)
    print("time taken in %s hrs" %((timer.time()-t0)/3600))
    rt.kill_phrqc_workers()
    #%%plot results
    plt.figure()
    plt.plot(time,AvgCa)
    plt.xlabel('Time [s]')
    plt.ylabel('Avg. Ca conc in aqeuous phase [mM]')
    plt.show()
    plt.figure()
    plt.imshow(rt.solid.portlandite.c)
    plt.show()
    