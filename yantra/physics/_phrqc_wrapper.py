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
#PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#This a wrapper to couple phreeqc with any transport code
#
#=======================================================================================

import numpy as np
import warnings
import sys
from copy import deepcopy  
import os
import multiprocessing as mp
from . import _multiprocessing_toolbox as _mt
from queue import Empty as QueueEmptyException

        
try:
   import IPhreeqcPy
except ImportError:
	warnings.warn("IPhreeqcPy not installed. Multicomponent reactive transport modules will not work.")

class Blank:  pass

class  ParallelPhrqcFunc:
    _args_name = ['output','dist_input_id']
    def __init__(self,*args):
        for name,val in zip(self._args_name,args):
            setattr(self,name,val)

def setvar(attr):
    """
    method to set property
    """
    def set_var(self,val):
        if self._vars[attr[1:]].type=='dict':
            try:
                oldval=getattr(self,attr)
                oldval.update(val)
                setattr(self,attr,oldval)
            except AttributeError:
                v = deepcopy(self._vars[attr[1:]].default)
                v.update(val)
                setattr(self,attr,v)
        elif self._vars[attr[1:]].type=='scalar':
            ones = np.ones(self.array_shape)
            val = val*ones
            if attr == '_poros':
                val[val>1]=1
                val[val<=0]=1
            setattr(self,attr,val.flatten(order='C'))
        else:
            setattr(self,attr,val)
    return set_var

def getvar(attr):
    """
    method to get property
    """
    def get_var(self):
        if self._vars[attr[1:]].type=='scalar':
            return getattr(self,attr).reshape(self.array_shape, order='C')
        else:
            return getattr(self,attr)
    return get_var

def delvar(attr):
    """
    method to delete property
    """
    def del_var(self):
        delattr(self,attr)
        return
    return del_var

def process_setvar(attr):
    """
    method to set property for PhrqcProcessWrapper
    """
    def set_var(self,val):
        self._run_cmd('set_var',[attr,val])
        self.get_output()
    return set_var
    
def process_getvar(attr):
    """
    method to get property for PhrqcProcessWrapper
    """
    def get_var(self):
        self._run_cmd('get_var',[attr])
        return self.get_output()
    return get_var
    
def process_delvar(attr):
    """
    method to delete property for PhrqcProcessWrapper
    """
    def del_var(self):
        self._run_cmd('del_var',[attr])
    return del_var

def parallel_setvar(attr):
    """
    method to set properties for ParallelPhrqc class
    """
    def set_var(self,val):
        if (attr in self.common_attrs) or (attr in self.additive_attrs):
            for i in range(self.nprocs):
                setattr(self.phrqc_workers[i],attr,val)
        elif attr in self.dist_attrs:
            val = self.splitter(val)
            for i,v in zip(list(range(self.nprocs)),val):
                setattr(self.phrqc_workers[i],attr,np.array(v))
        return 
    return set_var

def parallel_getvar(attr):
    """
    method to get properties for ParallelPhrqc class
    """
    def get_var(self):
        out = [0]*self.nprocs
        for i in range(self.nprocs):
            out[i]=getattr(self.phrqc_workers[i],attr)
        if attr in self.common_attrs:
            return out[0]
        elif attr in self.additive_attrs:
            return np.sum(np.array(out))
        elif (attr in self.dist_attrs):
            return self.merger(out)
    return get_var
    
def parallel_delvar(attr):
    """
    method to delete properties for ParallelPhrqc class
    """
    def del_var(self):
        for i in range(self.nprocs):
            delattr(self.phrqc_workers[i],attr)
        return
    return
            
def call_phrqc_method(name):
    def run_cmd(self,*args):
        return self._run_cmd(name,args)
    return run_cmd

def call_parallel_phrqc_method(name):
    """
    function to call methods mapped in ParallelPhrqc class
    """
    def call_func(self,*args):
        args= list(args)
        minfo = self.mapped_func[name]
        for i in range(len(args)):
            if i in minfo.dist_input_id:
                args[i]=self.splitter(args[i])
        for i in range(self.nprocs):
            loc_args = []
            for j in range(len(args)):
                if j in minfo.dist_input_id:
                    loc_args.append(args[j][i])
                else:
                    loc_args.append(args[j])
            getattr(self.phrqc_workers[i],name)(*loc_args)
        out = [0]*self.nprocs
        for i in range(self.nprocs):
            out[i] =self.phrqc_workers[i].get_output()
        if minfo.output == 'same':
            return out[0]
        elif minfo.output == 'merge':
            return self.merger(out)         
    return call_func
    
def _get_nprocs():
    try:
        return int(os.environ['OMP_NUM_THREADS'])
    except KeyError:
        return mp.cpu_count()

def dummy_domain(d):
    dummy =Blank()
    if d.d==2:
        dummy.nx,dummy.ny,dummy.d = d.nx,d.ny,d.d 
    elif d.d ==3:        
        dummy.nx,dummy.ny,dummy.nz,dummy.d = d.nx,d.ny,d.nz,d.d         
    return dummy
    
def phrqc(domain,domain_params,bc_params,solver_params):
    nprocs = _get_nprocs()
    if nprocs==1:
        return Phrqc(domain,domain_params,bc_params,solver_params)
    else:
        return ParallelPhrqc(domain,domain_params,bc_params,solver_params)

class Variable:
    def __init__(self,**kwargs):
        for k,v in kwargs.items():
            setattr(self,k,v)
            
class PhrqcMeta(type):
    def __init__(cls, name, bases, dct):
        """
        meta class of the LB created to initalize physical variable as property
        """
        for k,v in cls._vars.items():
            setattr(cls, k,property(fset=setvar("_"+k), fget=getvar("_"+k),
                             fdel=delvar("_"+k)))

class ParallelPhrqcMeta(type):
    """
    meta class for ParallelPhrqc class
    """
    def __init__(cls,name,bases,dct):
        for l in (cls.common_attrs,cls.dist_attrs,cls.additive_attrs):
            for v in l: 
                setattr(cls,v,property(fset=parallel_setvar(v),fget=parallel_getvar(v),
                                       fdel=parallel_delvar(v)))
            for m in cls.mapped_func.keys():
                setattr(cls,m,call_parallel_phrqc_method(m))       

class PhrqcProcessMeta(type):
    def __init__(cls, name, bases, dct):
        """
        meta class for PhrqcProcessWrapper class
        """
        for v in cls.mapped_attrs:
            setattr(cls, v,property(fset=process_setvar(v), fget=process_getvar(v),
                             fdel=process_delvar(v)))
        for m in cls.mapped_func:
            setattr(cls,m,call_phrqc_method(m))       

class ParallelPhrqc(object, metaclass=ParallelPhrqcMeta):
    """
    parallel wrapper that couples phreeqc with transport codes
    """
    common_attrs=['phrqc_input','selected_output_str','components','active_components',
                 'database', 'phrqc_input_file','boundary_solution_labels',
                 'eq_names','ss_names','kin_names','tracer_components','O_norm','H_norm',
                 'extra_selected_out_punch_head','extra_selected_out_punch_statements','time',
                 'dt','iters','phrqc_smart_run_tol','phrqc_flags','prev_selected_output']
    dist_attrs=['component_conc','solid_phase_conc','dphases','solution_labels','poros']
    additive_attrs=['nactive']
    mapped_func={'boundary_conditions':ParallelPhrqcFunc('same',[]),
                 'initial_conditions':ParallelPhrqcFunc('merge',[]),
                 'modify_solution':ParallelPhrqcFunc('merge',[0,2]),
                 'selected_output':ParallelPhrqcFunc('merge',[]),
                 'modify_solid_phases':ParallelPhrqcFunc('common',[0,1,2,3,4,5])
                 }
                 
    def __init__(self, domain,domain_params,bc_params,solver_params):
        mp_split_type=solver_params.get('mp_split_type','x_split')
        mp_chunks=solver_params.get('mp_chunks',1)
        self._splitter = getattr(_mt,mp_split_type)
        self._merger =  getattr(_mt,_mt.get_merge_type(mp_split_type))
        self.nprocs = mp_chunks* _get_nprocs()
        self.phrqc_workers=[]
        if domain.d == 2:
            array_shape= (domain.ny,domain.nx)
        elif domain.d == 3:
            array_shape = (domain.nz,domain.ny,domain.nx)
        slabels=domain_params.get('solution_labels',0)
        slabels= np.ones(array_shape)*slabels
        slabels=self.splitter(slabels)
        for p in range(self.nprocs):
            do = Blank()
            do.d = domain.d
            if domain.d == 2:
                do.ny,do.nx = np.shape(slabels[p])
            elif domain.d == 3:
                do.nz,do.ny,do.nx = np.shape(slabels[p])
            dp=deepcopy(domain_params)
            dp['solution_labels']=slabels[p]
            self.phrqc_workers.append(PhrqcProcessWrapper(p,do,dp,bc_params,solver_params))
    
    def kill(self):
        for worker in self.phrqc_workers:
            worker.kill()
    
    def splitter(self,val):
        if type(val).__name__ == 'dict':
            all_arrays=np.all([type(v).__name__=='ndarray' for v in list(val.values())])
            all_lists=np.all([type(v).__name__=='list' for v in list(val.values())])
            if all_arrays or all_lists:
                arranged_out=[]                 
                for p in range(self.nprocs):
                    arranged_out.append({})
                for k,v in val.items():
                    out= self.splitter(v)
                    for p in range(self.nprocs):                 
                        arranged_out[p][k]=out[p]
                return arranged_out
            else:
                return {} 
        elif type(val).__name__ == 'list':
#            all_arrays=np.all([type(v).__name__=='ndarray' for v in val])
#            all_lists=np.all([type(v).__name__=='list' for v in val])
#            if all_arrays:
            arranged_out = []
            for p in range(self.nprocs):
                arranged_out.append([])
            for v in val:
                out = self.splitter(v)
                for i in range(self.nprocs):
                    arranged_out[i].append(out[i])
            return arranged_out
        elif type(val).__name__=='ndarray':
            return self._splitter(val,self.nprocs)
        else:
            []
            
    def merger(self,val):
        a_condition=all([type(v).__name__=='dict' for v in val])
        if a_condition:
            out={}
            for v in val:
                for k in v.keys():
                    if k in out:
                        out[k].append(v[k])
                    else:
                        out[k]=[]
                        out[k].append(v[k])
            for k in list(out.keys()):
                out[k]=self._merger(out[k])
            return out
        else:
            return self._merger(val)
                    
class PhrqcProcessWrapper(object, metaclass=PhrqcProcessMeta):
    mapped_attrs=['phrqc_input','selected_output_str','components','active_components',
                 'component_conc','solid_phase_conc','dphases','nactive','database',
                 'phrqc_input_file','startcell','stopcell','solution_labels','poros',
                 'boundary_solution_labels','eq_names','ss_names','kin_names','tracer_components',
                 'O_norm','H_norm','extra_selected_out_punch_head',
                 'extra_selected_out_punch_statements','time','dt','iters','phrqc_smart_run_tol',
                 'phrqc_flags','prev_selected_output']
    mapped_func=['boundary_conditions','initial_conditions','modify_solution','active_nodes',
                 'modify_eq','modify_kin','modify_ss','tomoles','selected_output','sort_phases',
                 'modify_solid_phases','flatten_dict','reshape_dict','ndarray_dict',
                 'list_dict','add_dict','set_var','get_var','del_var']
                   
    def __init__(self,pid,domain,domain_params,bc_params,solver_params):
        self.pid = pid
        self.mp_timeout = solver_params.get('mp_timeout',30)
        self.queues = {'in':mp.Queue(),'out':mp.Queue(),'cmd':mp.Queue()}
        self.process = mp.Process(target = self._task, args = (self.queues,
                                                            self.mp_timeout))
        self.process.start()
        self.killed=False
        self._init_phrqc(dummy_domain(domain),domain_params,bc_params,solver_params)

    def kill(self):
        self.queues['cmd'].put('kill')
        self.process.join()
        
    def _run_cmd(self,cmd,args=[]):
        mapped_func=self.mapped_func
        if not self.killed:
            if cmd in mapped_func:
                self.queues['cmd'].put(cmd)
                self.queues['in'].put(args)
#                return self.get_output()
            else:
                raise ValueError('only functions can be ran with this command')
        else:
            raise ValueError("process have been killed-calculations can't be performed anymore")
            
    def _task(self, queues, timeout):
        """
        phreeqc worker is initialized in this function for the process
        """
        mapped_func,mapped_attrs=self.mapped_func,self.mapped_attrs
        q = queues
        cmd = q['cmd'].get()
        while cmd != 'kill':            
            if cmd == 'init':
                phrqc = Phrqc(*q['in'].get(timeout=timeout))
                q['out'].put(1)
            else:
                if cmd in mapped_func:
                    q['out'].put(getattr(phrqc,cmd)(*q['in'].get(timeout=timeout)))
                elif cmd in mapped_attrs:
                    q['out'].put(getattr(phrqc,cmd))
            try:
                # timeout is necessary to prevent the process from clogging
                # up the CPU
                cmd = q['cmd'].get()
                if cmd=='kill':
                    self.killed = True
            except:
                pass
#            except QueueEmptyException:
                # If a process hasn't been commanded to do anything within timeout,
                #we can assume that the parent process is finished.
                # This process will end itself to preserve memory
#                cmd = 'kill'
#                self.killed=True
    
    def get_output(self):
       out= self.queues['out'].get()
       try:
            while True:
                self.queues['out'].get_nowait()
       except QueueEmptyException:
            pass
       return out
   
    def _init_phrqc(self,domain,domain_params,bc_params,solver_params):
        self.queues['cmd'].put('init')
        self.queues['in'].put((domain,domain_params,bc_params,solver_params))
        return self.get_output()
        
        

class Phrqc(object, metaclass=PhrqcMeta):

    """
    wrapper to couple phreeqc with transport codes
    """
    _vars={'database':Variable(default='phreeqc.dat',type='str',loc='domain_params',
                               doc='path for database'),
           'phrqc_input_file':Variable(default='',type='str',loc='domain_params',
                                doc='path for phreeqc input file'),
           'startcell':Variable(default=0,type='int',loc='solver_params',
                                doc='number of start cell'),
           'stopcell':Variable(default=0,type='int',loc='solver_params',
                               doc='number of end cell'),
           'solution_labels':Variable(default=0,type='scalar',loc='domain_params',
                                      doc='array of solution lable cells are assigned'),
           'poros':Variable(default=1,type='scalar',loc='domain_params',
                                      doc='water filled porosity'),
           'boundary_solution_labels':Variable(default={},type='dict',
                                                loc='boundary_params',doc='solution labels for boundary'),
           'eq_names':Variable(default=[],type='list',loc='solver_params',
                                    doc='list of names of equilibrium phases'),
           'ss_names':Variable(default={},type='dict',loc='solver_params',
            doc='dictonary listing solid solution name as keys and each keys corresponds to link of components in solid solution'),
           'kin_names':Variable(default=[],type='list',loc='solver_params',
                                    doc='list of names kinetic phases'),
           'tracer_components':Variable(default=[],type='list',loc='solver_params',
                                        doc='list of tracer components'),
           'O_norm':Variable(default=55.506216797268586,type='float',
                             doc='number of moles of O in water'),
           'H_norm':Variable(default=111.01243359454628,type='float',loc='solver_params',
                             doc='number of moles of H in water'),
           'extra_selected_out_punch_head':Variable(default=[],type='list',loc='solver_params',
                                                    doc='headings for extra punch statements in selected output'),
           'extra_selected_out_punch_statements':Variable(default=[],type='list',loc='solver_params',
                                                          doc = 'punch statements for selected output'),
           'time':Variable(default=0,type='float',loc='solver_params',doc='time till simulation'),
           'dt':Variable(default=0,type='float',loc='solver_params',doc='timestep'),
           'iters':Variable(default=0,type='int',loc='solver_params',doc='number of iterations'),
           'nactive':Variable(default=0,type='int',loc='solver_params',doc='number of active phreeqc cells'),
           'phrqc_smart_run_tol':Variable(default=1e-8,type='float',loc='solver_params',doc='tolerance for smart run'),
           'phrqc_flags':Variable(default={'only_interface':False,
                                           'only_fluid':False,
                                           'smart_run':False,
                                           'include_OH':True,
                                           'update_output':True,
                                           'update_poros':False,
                                            },type='dict',loc='solver_params', 
                                            skip_initalization=True,
                                            doc='flags for phrqc wrapper'),
         'prev_selected_output':Variable(default={},type='dict',loc='domain_params',
                                    doc='private variable for phreeqc selected output for previous timestep'),
         'ss':Variable(default={},type='dict',loc='domain_params',
                                    doc='source sink term for transport step'),
        }

    
    def __init__(self, domain,domain_params,bc_params,solver_params):
        """
        Initialize phreeqc model

        Arguments

        """
        self.nprocs = 1
        self._selected_output ={}
        if domain.d == 2:
            self.array_shape= (domain.ny,domain.nx)
        elif domain.d == 3:
            self.array_shape = (domain.nz,domain.ny,domain.nx)
        else:
            raise ValueError('dimension of the domain not understood')
        inputs = deepcopy(domain_params)
        inputs.update(solver_params)
        if 'solution_labels' in bc_params:
            inputs['boundary_solution_labels']= bc_params['solution_labels']
        for k,v in self._vars.items():
            setattr(self,k,inputs.get(k,v.default))
        self.ncells= len(self._solution_labels)
        try:
            assert self.ncells == (self.stopcell-self.startcell)+1
        except AssertionError:
            self.stopcell = self.ncells
            self.startcell = 1
        self.IPhreeqc = IPhreeqcPy.IPhreeqc()
        self.IPhreeqc.LoadDatabase(self.database)
        self.IPhreeqc.RunFile(self.phrqc_input_file)
        self.IPhreeqc.RunString(self.phrqc_input)
        self.phrqc_flags['update_output']=True
        self.O_norm = self.selected_output()['H2O'][0]
        self.H_norm = 2.0 * self.selected_output()['H2O'][0]

    @property
    def phrqc_input(self):
        """
        input for phreeqc includes selected output string
        """
        with open(self.phrqc_input_file,'r') as f:
            string = f.read()
        string +='\n'+self.selected_output_str
        string +='\nend'
        return string
    
    @property
    def selected_output_str(self):
        """
        Generates selected output string
        """
        selected_output= []
        selected_output.append('SELECTED_OUTPUT')
        selected_output.append('\t-reset false')
        selected_output.append('\t-time false')
        selected_output.append('\t-high_precision true')
        selected_output.append('\t-solution true')
        selected_output.append('\t-pH true')
        selected_output.append('\t-pe false')
        selected_output.append('\t-charge_balance false')
        selected_output.append('\t-alkalinity true')
        selected_output.append('\t-ionic_strength true')
        selected_output.append('\t-percent_error false')
        selected_output.append('USER_PUNCH')
        user_punch_head =[]
        user_punch_str = []
        counter = 0
        user_punch_head.append("H2O")
        counter+=10
        user_punch_str.append(str(counter)+'\tpunch'+'\tmol("H2O")')
        user_punch_head.append("poros")
        counter+=10
        user_punch_str.append(str(counter)+'\tpunch'+'\ttot("water")')
        for name in self.components:
            user_punch_head.append(name)
            counter+=10
            user_punch_str.append(str(counter)+'\tpunch'+'\ttot("%s")'%name)
        for name in self.eq_names:
            user_punch_head.append(name)
            counter+=10
            user_punch_str.append(str(counter)+'\tpunch'+'\tequi("%s")'%name)
            user_punch_head.append('SI_'+name)
            counter+=10
            user_punch_str.append(str(counter)+'\tpunch'+'\tSI("%s")'%name)
        for name in self.kin_names:
            user_punch_head.append(name)
            counter+=10
            user_punch_str.append(str(counter)+'\tpunch'+'\tKIN("%s")'%name)
            user_punch_head.append('SI_'+name)
            counter+=10
            user_punch_str.append(str(counter)+'\tpunch'+'\tSI("%s")'%name)
        ss_component_names = []
        for comp in list(self.ss_names.values()):
            ss_component_names.extend(comp)
        for name in ss_component_names:
            user_punch_head.append(name)
            counter+=10
            user_punch_str.append(str(counter)+'\tpunch'+'\ts_s("%s")'%name)
            user_punch_head.append('SI_'+name)
            counter+=10
            user_punch_str.append(str(counter)+'\tpunch'+'\tSI("%s")'%name)
        user_punch_head.extend(self.extra_selected_out_punch_head)
        counter+=10
        extra_statements=str(counter)+'\tpunch'+'\t'.join(self.extra_selected_out_punch_statements)
        user_punch_str.append(extra_statements)
        selected_output.append('\t-headings %s'%(' '.join(user_punch_head)))
        selected_output.append('\t-start')
        selected_output.append('\n\t\t'.join(user_punch_str))
        selected_output.append('\t-end')
        return '\n'.join(selected_output)
    
    @property 
    def components(self):
        """
        list of components in phreeqc
        """
        try:
            getattr(self,'_components')
        except AttributeError:
            components= self.IPhreeqc.GetComponentList()
            if 'H' not in components:
                components.append('H')
            if 'O' not in components:
                components.append('O')
            if not self.phrqc_flags['include_OH']:
                components.remove('H')
                components.remove('O')
            self._components=components
        return self._components
    
    @property
    def active_components(self):
        """
        lists only non tracer components
        """
        return [name for name in self.components
                                if name not in self.tracer_components]
    
    def  boundary_conditions(self):
        """
        gets boundary condition for transport model
        """
        if self.time == 0:
            results = self.selected_output()
            self._bc_solutions = {}
            for bc,label in self.boundary_solution_labels.items():
                self._bc_solutions[bc]={}
                idx = results['soln'].index(label)
                for  name in self.components:
                    self._bc_solutions[bc][name]=results[name][idx]
        return self._bc_solutions
   
    def initial_conditions(self):
        """
        gets initial conditions for transport model
        """
        if self.time == 0:
            copycellstr = []
            for (cell, label) in zip(list(range(self.startcell,self.stopcell+1,1)),
                 self._solution_labels):
                copycellstr.append('copy cell %i %i' % (label, cell))
            copycellstr = '\n'.join(copycellstr)
            runcellstr=['RUN_CELLS']
            runcellstr.append('\t-cells %i-%i'%(self.startcell,self.stopcell))
            runcellstr.append('\t-start_time 0')
            runcellstr.append('\t-time_step 0')
            runcellstr = '\n'.join(runcellstr)
            self.IPhreeqc.AccumulateLine(copycellstr)
            self.IPhreeqc.AccumulateLine('end')
            self.IPhreeqc.AccumulateLine(runcellstr)
            self.IPhreeqc.AccumulateLine('end')
            self.IPhreeqc.RunAccumulated()
            self.phrqc_flags['update_output'] = True
            self.selected_output(merge_with_previous=False,toarray=True,reshape=True)
            self.poros = self.selected_output()['poros']
            self._initial_conditions=deepcopy(self.component_conc)
        return self._initial_conditions
                
    def modify_solution(self,c_trans,dt,nodetype):
        """
        inconc/prevconc flatten dictonary of conc subset.
        nodetype/solid_phase_qty:matrix or flattened either is possible
        to be run by current iphreeqc solver.

        input:
        -----
        inconc = dictionary  containing each component as key and ncells \
        values of conc for each key
        current_time=current time in sec
        nodetype= node type matrix

        output:
        -------
        updates currentselecout
        """
        self.dt = dt
        active_nodes = self.active_nodes(c_trans,nodetype)
        moles = self.flatten_dict(self.to_moles(c_trans))
        modify_str=[]
        runcells=[]
        runcell_str=[]
        for i,cell in enumerate(range(self.startcell,self.stopcell+1,1)):
            if active_nodes[i]:
                runcells.append(str(cell))
                modify_str.append('SOLUTION_MODIFY %i' % cell)
                if 'H' in moles:
                    c = moles['H'][i]
                    modify_str.append('\t%s\t%.20e' % ('total_h', c))
                if 'O' in moles:
                    c = moles['O'][i]
                    modify_str.append('\t%s\t%.20e' % ('total_o', c))
                modify_str.append('\t-totals')
                for name,val in moles.items():
                    if (name != 'H') and (name!='O'):
                        c = val[i]
                        modify_str.append('\t\t%s\t%.20e' % (name, c))
        modify_str.append('end')
        modify_str ='\n'.join(modify_str)
        runcell_str.append('RUN_CELLS')
        runcell_str.append('\t-time_step %s'%self.dt)
        runcell_str.append('\t-cells %s'%'\n\t\t'.join(runcells))
        runcell_str.append('end')
        runcell_str='\n'.join(runcell_str)
        self.IPhreeqc.AccumulateLine(modify_str)
        self.IPhreeqc.AccumulateLine(runcell_str)
        self.IPhreeqc.RunAccumulated()
        self.phrqc_flags['update_output'] = True
        self.selected_output(merge_with_previous=True,toarray=True,reshape=True)
        if self.phrqc_flags['update_poros']: self.poros = self.selected_output()['poros']
        c_current = self.component_conc
        ss={}
        for name in self.components:
            ss[name] = (c_current[name]-c_trans[name])/self.dt
            ss[name] = ss[name] *(active_nodes.reshape(self.array_shape)>0)
            ss[name] *= self.poros
        self.time+=self.dt
        self.iters+=1
        return ss #returns source-sink in physical units

    def active_nodes(self,c,nodetype):
        """
        gives active nodes for phreeqc
        """
        active =np.ones(self.array_shape).flatten()
        prev_c = self.component_conc
        smart_active = np.zeros(self.array_shape)
        inactive=np.ones(self.array_shape).flatten()
        tot_solid_phase_conc = self.add_dict(self.flatten_dict(self.solid_phase_conc))
        self.tot_solid_phase_conc=tot_solid_phase_conc
        inactive-=1*np.logical_and(tot_solid_phase_conc > 1e-14, nodetype.flatten(order='c')<=0)
        if self.phrqc_flags['only_interface']:
            inactive -=1*(nodetype.flatten()==0) #activate interface nodes
        elif self.phrqc_flags['only_fluid']:
            inactive-=1*(nodetype.flatten()<=0) #activate liquid nodes
        elif self.iters >1 and self.phrqc_flags['smart_run']:
            for name in self.active_components:
                diff =np.abs( c[name]*(c[name]-prev_c[name])/(c[name]+1e-30)**2)
                smart_active +=1*(diff<self.phrqc_smart_run_tol)
            inactive-=1*(smart_active.flatten(order='c') <=0)
        else:
            inactive = np.zeros(self.array_shape).flatten()
        active -= 1*(inactive>0)
        self.nactive = np.sum(active)
        return active
            
    @property    
    def component_conc(self):
        """
        concentration of aqueous components 
        """
        selected_output = self.selected_output()
        output={}
        for name in self.components:
            output[name]=selected_output[name]
        return output
    
    @property
    def solid_phase_conc(self):
        """
        concentration of solid phases
        """
        selected_output = self.selected_output()
        phases={}
        for name in self.eq_names:
            phases[name]= selected_output[name]
        for name in self.kin_names:
            phases[name]= selected_output[name]
        for components in self.ss_names.values():
            for component in components:
                phases[component]= selected_output[component]
        return phases
                      
    def modify_eq(self, phaseqty, si, dissolve_only, precipitate_only, skip_phase):
        """
        modifies the phaseqty in phreeqc
        
        Parameters
        ----------
        phaseqty: dict
            dictionary containing an ndarray of new quantities of equilibrium phases
        """
        phaseqty = self.flatten_dict(phaseqty)
        si = self.flatten_dict(si) 
        dissolve_only = self.flatten_dict(dissolve_only)
        precipitate_only = self.flatten_dict(precipitate_only)
        skip_phase = self.flatten_dict(skip_phase)
        modifystr = []
        for cell in range(self.startcell,self.stopcell+1,1):
            modifystr.append("EQUILIBRIUM_PHASES_MODIFY %d" % cell)
            for key in list(phaseqty.keys()):
                if key in skip_phase:
                    cdn =1- (skip_phase[key][cell-1]<=0)
                else:
                    cdn = 1
                if cdn:                        
                    modifystr.append("\t -component %s" %(key))
                    c =phaseqty[key][cell-1]
                    if c < 0: c=1e-30
                    modifystr.append("\t\t%s\t%.20e" %('-moles', c))
                    if key in si:
                        modifystr.append("\t\t%s\t%.20e" %('-si', si[key][cell-1]))
                    if key in dissolve_only:
                        modifystr.append("\t\t%s\t%i" %('-dissolve_only', dissolve_only[key][cell-1]))
                    if key in precipitate_only:
                        modifystr.append("\t\t%s\t%i" %('-precipitate_only', precipitate_only[key][cell-1]))
        modifystr.append("end")       
        modifystr ='\n'.join(modifystr)
        self.IPhreeqc.AccumulateLine(modifystr)

    def modify_kin(self, kinqty, d_params):
        """
        modifies the kinetic phases in phreeqc
        
        Parameters
        ----------
        kinqty: dict
            dictionary containing an ndarray of new quantities of kinetic phases
        """
        kinqty = self.flatten_dict(kinqty)
        for k,v in d_params.items():
            for i in range(len(v)):
                d_params[k]= [val.flatten(order='c') for val in v]                    
        modifystr = []
        for cell in range(self.startcell,self.stopcell+1,1):
            modifystr.append("KINETICS_MODIFY %d" % (cell))
            for key in list(kinqty.keys()):
                modifystr.append("\t -component %s" %key)
                m = kinqty[key][cell-1]
                if m < 0: m=1e-30                    
                modifystr.append("\t\t%s\t%.20e" %('-m', m))
                modifystr.append("\t\t%s\t%.20e" %('-m0',m))
                if key in d_params:
                    dpar = [d_params[key][i][cell-1] for i in range(len(d_params[key]))]
                    str_d_params = ' '.join(str(i) for i in dpar)
                    modifystr.append("\t\t%s\n\t\t\t%s" %('-d_params', str_d_params))
        modifystr.append("end")       
        modifystr ='\n'.join(modifystr)
        self.IPhreeqc.AccumulateLine(modifystr)
             
    def modify_ss(self,ssqty):
        """
        modify solid solution in phreeqc
        
        Parameters
        ----------
        ssqty: dict
            dictionary containing an ndarray of new quantities of kinetic phases
        """
        ssqty = self.flatten_dict(ssqty)
        modifystr = []
        for cell in range(self.startcell,self.stopcell+1,1):
            modifystr.append("SOLID_SOLUTIONS_MODIFY %d" % (cell))
            for key in list(ssqty.keys()):
                if key in self.ss_names:
                    modifystr.append("\t-solid_solution %s"%key)
                    for comp in list(ssqty[key].keys()):
                        modifystr.append("\t-component %s"%comp)
                        m = ssqty[key][comp][cell-1]
                        if m <0: m = 1e-30 
                        modifystr.append("\t\t%s\t%.20e" %('-moles',m))
        modifystr.append("end")       
        modifystr ='\n'.join(modifystr)
        self.IPhreeqc.AccumulateLine(modifystr)
        
    def to_moles(self, c):
        """
        converts concentration to moles which can be used in modify_solution method
        
        Parameters
        ----------
        c: dict
            dictonary containing ndarray of cocentrations for different components 
        """
        moles = self.flatten_dict(c)
        for name in self.components:
            moles[name][np.isnan(moles[name])]=1e-30
            moles[name][moles[name]<0]=1e-30
            if name == 'H':
                moles[name] += self.H_norm
                moles[name] *= self._poros
            elif name == 'O':
                moles[name] += self.O_norm
                moles[name] *= self._poros
            else:
                pass
                moles[name] *= self._poros
        return moles

    def selected_output(self,merge_with_previous=False,toarray=False,reshape=False):
        """
        returns selected output. The new selected output is generated if phrqc_flags['update_output'] \
        is true else returns the currently saved output
        
        Parameters
        ----------
        merge_with_previous: boolean
            if true merges the current selected output with previously saved output
        
        toarray: boolean
            if true converts output into array
        
        reshape: boolean
            if true resphases arrays to the array_shape
        """
        if self.phrqc_flags['update_output']:
            self.prev_selected_output = deepcopy(self._selected_output)
            output=self.IPhreeqc.GetSelectedOutputArray()
            selected_output={}
            if len(output) > 0:
                header = output[0]
                for head in header:
                     selected_output[head] = []
                for row in output[1:]:
                    for (i, head) in enumerate(header):
                        results = row[i]
                        if head == 'H':
                            results -= self.H_norm
                        elif head == 'O':
                            results -= self.O_norm
                        selected_output[head].append(results)
            if merge_with_previous: 
                merged_output = self.list_dict(self.flatten_dict(self.prev_selected_output))
                for name,val in selected_output.items():
                    for i,cell in enumerate(selected_output['soln']): 
                        merged_output[name][cell-1]=val[i]             
                selected_output = merged_output
#           convert to ndarray
            if toarray:
                selected_output = self.ndarray_dict(selected_output)
                if reshape:
                    selected_output = self.reshape_dict(selected_output,self.array_shape)
            self._selected_output = selected_output
            self.phrqc_flags['update_output'] = False
        return self._selected_output
    
    @property
    def dphases(self):
        """
        difference of phases in a given timestep after reaction
        """
        output={}
        if self.time >0:
            current_selected_out = self.selected_output()
            prev_selected_out = self.prev_selected_output
            for name in self.eq_names:
                output[name] = current_selected_out[name] - prev_selected_out[name]
            for name in self.kin_names:
                output[name] = current_selected_out[name] - prev_selected_out[name]
            for ss in list(self.ss_names.keys()):
                for name in self.ss_names[ss]:
                    output[name] = current_selected_out[name] - prev_selected_out[name]
        return output
        
    def sort_phases(self,phaseqty):
        """
        sorts phases in equillibrium phases, solid solutions and kinetic phases
        
        Parameters
        ----------
        phaseqty: dict
            dictonary of ndarray giving quantities of phases
            
        Returns
        -------
        (dict, dict, dict)
            dictonary of eqillibrium phases, solid solutions and kinetic phases
        """
        eqphases = {}
        for name in self.eq_names:
            if name in phaseqty:
                eqphases[name] = phaseqty[name]
        ssphases = {}
        for name,components in self.ss_names.items():
            ssphases[name] ={}
            for component in components: 
                if component in phaseqty:
                    ssphases[name][component] = phaseqty[component]
        kinphases = {}
        for name in self.kin_names:
            if name in phaseqty:
                kinphases[name] = phaseqty[name] 
        return eqphases,ssphases,kinphases
                
    def modify_solid_phases(self,phaseqty,modifyeq_si,modifyeq_dissolve_only,modifyeq_precipitate_only,modifyeq_skip_phase,modifykin_d_params):
        """
        modifies solid phases to the qunatities given as input in phaseqty
        
        Parameters
        ----------
        phaseqty: dict
            dictonary of ndarray giving quantities of solid phases
        """
        #first sort the phases
        for name,val in phaseqty.items():
            self._selected_output[name] = val
        eqphases, ssphases, kinphases = self.sort_phases(phaseqty)
        if len(eqphases)>0:
            self.modify_eq(eqphases,modifyeq_si,modifyeq_dissolve_only,modifyeq_precipitate_only,modifyeq_skip_phase)
        if len(ssphases)>0:
            self.modify_ss(ssphases)            
        if len(kinphases)>0:
            self.modify_kin(kinphases,modifykin_d_params)     
        runcellstr=['RUN_CELLS']
        runcellstr.append('\t-cells %i-%i'%(self.startcell,self.stopcell))
        runcellstr = '\n'.join(runcellstr)
        self.IPhreeqc.AccumulateLine(runcellstr)
        self.IPhreeqc.AccumulateLine('end')
        self.IPhreeqc.RunAccumulated()
        self.phrqc_flags['update_output'] = True
        self.selected_output(merge_with_previous=True,toarray=True,reshape=True)
    
    def save_dict(self):
        out= {}
        for k in self._vars.keys():
            out[k] = getattr(self,k)
        out['dumpstr'] = self.IPhreeqc.GetDumpString()
         
        return out

    def set_var(self,attr,val):
        setattr(self,attr,val)
    
    def get_var(self,attr):
        return getattr(self,attr)
    
    def del_var(self,attr):
        return delattr(self,attr)

    def flatten_dict(self,d):
        """
        flattens ndarray in the dictionary
        
        Parameters
        ----------
        d: dict
            input dictionary
        """
        new_d = deepcopy(d)
        for k,v in new_d.items():
            if type(v).__name__ =='ndarray':
                new_d[k]=v.flatten(order='C')
            elif type(v).__name__=='dict':
                new_d[k]=self.flatten_dict(v)
        return new_d
    
    def kill(self):
        pass
    
    def reshape_dict(self,d,shape):
        """
        reshapes ndarray in the dictionary
        
        Parameters
        ----------
        d: dict
            input dictionary
        """
        new_d = d
        for k,v in new_d.items():
            if type(v).__name__ =='ndarray':
                new_d[k]=v.reshape(shape, order='C')
            elif type(v).__name__=='dict':
                new_d[k]=self.reshape_dict(v)
        return new_d
 
    @staticmethod
    def ndarray_dict(d):
        """
        convert list  dictionary into ndarray dictionary
        
        Parameters
        ----------
        d: dict
            input dictionary
        """
        new_d = d
        for k,v in new_d.items():
            new_d[k]=np.array(v)
        return new_d

    @staticmethod
    def list_dict(d):
        """
        convert ndarray  dictionary into list dictionary
        
        Parameters
        ----------
        d: dict
            input dictionary
        """
        new_d = d
        for k,v in new_d.items():
            new_d[k]=list(v)
        return new_d
    
    @staticmethod
    def add_dict(d):
        tot=0
        for v in d.values():
            tot += v
        return tot
