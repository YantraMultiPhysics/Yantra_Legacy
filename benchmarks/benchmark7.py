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
Benchmark 7: Phreeqc manual Example 11 is replicated in this benchmark. The results  \
of this benchmark are published in Patel et al (2013) coupled Problems conference paper.
"""
#required as phreeqc parallelization is done in yantra using python multiprocessing toolbox
#which requires this  line
if __name__ =='__main__':
    #%%import modules
    import sys
    import os
    PARENT = '..'
    sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
    import yantra
    from yantra import timer
    import numpy as np
    import matplotlib.pylab as plt
    #%% function to read comsol results
    def read_comsol_results():
        res=np.loadtxt('ext_results/benchmark7.csv',delimiter=',',skiprows=1)
        with open('ext_results/benchmark7.csv','r') as f:
            head=f.readline()
            head=head.split(',')
            head = [name.rstrip() for name in head]
        comsol={}
        for i,name in enumerate(head):
            comsol[name]=res[:,i]
        return comsol
    #%%set domain
    domain = yantra.Domain2D((0,0),(0.08,0.04),0.08/40.,grid_type='midway')
    #%%set reactive transport model
    #domain params
    ux=0.00027777/100
    domain_params={}
    domain_params['D']=(0.2/100)*(0.00027777/100)
    domain_params['u']=[ux,0]
    domain_params['database']='phreeqc.dat'
    domain_params['phrqc_input_file']='benchmark7.phrq'
    domain_params['solution_labels']=100002
    domain_params['tracer_components']=['Cl','N']
    #bc params
    bc_params={}
    bc_params['solution_labels']={'left':100001}
    bc_params['top']=['flux',0]
    bc_params['bottom']=['flux',0]
    bc_params['left']=['c',0]
    bc_params['right']=['open',0]
    #solver parameters
    solver_params={}
    solver_params['collision_model']='srt'
    solver_params['phrqc_flags']={}
    solver_params['phrqc_flags']['smart_run']=False
    solver_params['phrqc_smart_run_tol']=1e-6
    solver_params['mp_split_type']='y_split'
    rt= yantra.PhrqcReactiveTransport('AdvectionDiffusion', domain,
                                      domain_params,bc_params,solver_params)
    #%%adjust bc for setting cauchy
    for name in rt.fluid.components:
        comp=getattr(rt.fluid,name)
        comp.bc={'left':['flux',comp.bc['left'][1]*ux]}
    #%%run model
    results={}
    time=[]
    clock=timer()
    for name in rt.fluid.components:
        if name not in ['H','O']:
            results[name]=[]
    while rt.time <=86400:
        rt.run(iters=rt.iters+5)
        print("time: %s"%rt.time)
        time.append(rt.time)
        for name in rt.fluid.components:
            if name not in ['H','O']:
                comp=getattr(rt.fluid,name)
                results[name].append(comp.c[int(comp.ny/2),-1])
    rt.kill_phrqc_workers()
    print('time taken: %s mins'%clock.time_lapsed)
    #%%plot_results
    plt.figure()
    #plot lb results
    for name,v in results.items():
        plt.plot(time,v,'--',label=name)
    #plot comsol_results
    comsol = read_comsol_results()
    for name,v in comsol.items():
        if name!='time':
            plt.plot(comsol['time'],v,'o',markevery=10,label=name)
    plt.legend()
    plt.xlabel('time (s)')
    plt.ylabel('c (mol/lit)')
    plt.show()
     
