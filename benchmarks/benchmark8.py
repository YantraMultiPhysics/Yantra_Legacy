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
Benchmark 8: Portlandite dissolution without geometry update. The results  \
of this benchmark are published in Patel et al (2014), Phy & Chem. of Earth.
"""
if __name__=='__main__':
    #%% import modules
    import sys,os
    PARENT = '..'
    sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
    import yantra
    import numpy as np
    import matplotlib.pylab as plt
    #%% read HP1 results
    def read_hp1_results():
        pass
    #%% set domain
    domain = yantra.Domain2D((0,0),(3e-2,1e-2),3e-2/80.)
    domain.nodetype[:,-1] = 1
    #%%set reactive transport model
    #domain params
    domain_params={}
    domain_params['D']=1e-9
    domain_params['database']='cemdata07.dat'
    slabels = 100002 *np.ones(domain.shape)
    slabels[:,-2] = 100003
    domain_params['phrqc_input_file']='benchmark8.phrq'
    domain_params['solution_labels']=slabels
    domain_params['eq_names']=['portlandite']
    #bc params
    bc_params={'solution_labels':{'left':100001},
               'left':['c',0],
               }
    #solver parameters
    solver_params={}
    solver_params['collision_model']='srt'
    solver_params['phrqc_flags']={}
    solver_params['phrqc_flags']['smart_run']=False
    #solver_params['phrqc_smart_run_tol']=1e-8
    solver_params['mp_split_type']='y_split'
    
    #generate model
    rt= yantra.PhrqcReactiveTransport('Diffusion', domain,
                                      domain_params,bc_params,solver_params)  
    
    #%%run model
    Ts=60*3600            
    results={}
    for name in ['time','portlandite','pH']:
        results[name]=[]
    while rt.time <=Ts:
        rt.advance()
        if (rt.iters%5==0):
            print('Current Time: %s' %rt.time)
            results['time'].append(rt.time/3600)
            results['pH'].append(rt.phrqc.selected_output()['pH'][int(domain.ny/2),-2])
            results['portlandite'].append(rt.phrqc.selected_output()['portlandite'][int(domain.ny/2),-2]*1000)
    rt.kill_phrqc_workers()
    #%%plot results  
    f, axarr = plt.subplots(1,2, sharex=False,sharey=False)
    axarr[0].set_title('Portlandite concentration')
    axarr[0].plot(results['time'], results['portlandite'],label='LB')
    axarr[0].set_xlabel('Time (hr)')
    axarr[0].set_ylabel('Concentration (mmol/lit)')
    axarr[1].set_title('pH')
    axarr[1].plot(results['time'], results['pH'],label='LB')
    axarr[1].set_ybound(12.3,12.5)
    axarr[1].set_xlabel('Time (hr)')
    axarr[1].set_ylabel('pH')
    plt.show()
            
            
            
            