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
#=======================================================================================

import numpy as np
import yantra
from yantra._base import Multicomponent 
from yantra._base import Fmethod
from ._phrqc_wrapper import  phrqc
#from copy import deepcopy 

class PhrqcReactiveTransport(object):
    """
    A Lattice Boltzmann method based reactive transport model where reactions are computed using \
    geochemical solver phreeqc
    """
    _signature = 'yantra.physics.PhrqcReactiveTransport.PhrqcReactiveTransport'
    def __init__(self,eqn,domain,domain_params,bc_params,solver_params):
        """
        A Lattice Boltzmann method based reactive transport model where reactions are computed using \
        geochemical solve phreeqc
        """
        self.auto_time_step = solver_params.get('auto_time_step',True)
        self.phrqc = phrqc(domain,domain_params,bc_params,solver_params)
        self.solid=Solid(domain,domain_params,solver_params)
        components = self.phrqc.components
        bc = self.phrqc.boundary_conditions() 
        init_c = self.phrqc.initial_conditions()
        if self.solid.nphases >0:
            phaseqty = self.solid.phaseqty_for_phrqc()
            self.phrqc.modify_solid_phases(phaseqty,self.modifyeq_si(),self.modifyeq_dissolve_only(),
                                           self.modifyeq_precipitate_only(),self.modifyeq_skip_phase(),
                                           self.modifykin_d_params())
        for name in components:
            if name not in domain_params:
                domain_params[name]={}
            domain_params[name]['c']=init_c[name]
        for name in bc:
            for comp in components:
                if comp not in bc_params:
                    bc_params[comp]={}
                bc_params[comp][name]=['c',bc[name][comp]]
        #tortousity model parameters
        self.app_tort_params=domain_params.get('app_tort_params',[1,0])
        #append parameters from solid in initialization 
        if  ('Multilevel' in eqn) and (self.solid.n_diffusive_phases>0):
            app_tort= self.app_tort()
            domain_params['app_tort'] = app_tort
            poros=self.solid.poros 
            poros[poros<=0]=1
            poros[poros>1]=1
            domain_params['poros']=poros
            #change phreeqc porosity
            setattr(self.phrqc,'poros',poros)        
        self.fluid=Multicomponent(eqn,components,domain,domain_params,
                                  bc_params,solver_params)

    def advance(self):
        """
        """
        #advance phrqc
        c=self.fluid.get_attr('c')
        ss=self.phrqc.modify_solution(c,self.dt,self.solid.nodetype)     
        self.fluid.set_attr('ss',ss) 
        #advance fluid
        self.fluid.call('advance') 
        #update solid properties
        poros = self.solid.poros
        phaseqty=self.solid.update(self.phrqc.dphases)
        dporosdt = (self.solid.poros-poros)/self.dt
        if len(phaseqty):
            self.phrqc.modify_solid_phases(phaseqty,self.modifyeq_si(),self.modifyeq_dissolve_only(),
                                           self.modifyeq_precipitate_only(),self.modifyeq_skip_phase(),
                                           self.modifykin_d_params())
            self.fluid.set_attr('nodetype',self.solid.nodetype,component_dict=False)
        if  ('Multilevel' in self.fluid.eqn) and (self.solid.n_diffusive_phases>0):
            poros=self.solid.poros
            poros[poros<=0]=1
            poros[poros>1]=1
            app_tort= self.app_tort()
            self.fluid.set_attr('dporosdt',dporosdt,component_dict=False) 
            self.fluid.call('update_transport_params',poros,
                            app_tort,self.auto_time_step)
            setattr(self.phrqc,'poros',poros)    

  
            
    def run(self,**kwargs):
        """
        runs an lb model for given time or number of iterations
        
        Parameters
        ----------
        time: float
            time till which simulation needs to carried out
        iters: int
            number of iterations till which simulations have to be carried out
            
        Note: Only keyword arguments allowed
        """
        if 'time' in kwargs:
            time = kwargs['time']
            while self.time < time:
                self.advance()
        elif 'iters' in kwargs:
            iters = kwargs['iters']
            while self.iters < iters:
                self.advance()
        else:
            raise ValueError("Either time or iters keywords can be specified")
    
    @property
    def dt(self):
        return getattr(self.fluid,self.fluid.components[0]).dt
    
    @property
    def time(self):
        return getattr(self.fluid,self.fluid.components[0]).time

    @property
    def iters(self):
        return getattr(self.fluid,self.fluid.components[0]).iters

    def app_tort(self):
        a,n = self.app_tort_params
        poros=self.solid.poros
        poros[poros<=0]=1
        poros[poros>1]=1
        return a* poros**n
    
    def modifyeq_si(self):
        return {}
    
    def modifyeq_skip_phase(self):
        return {}
    
    def modifykin_d_params(self):
        return {}
    
    def modifyeq_precipitate_only(self):
        return {}
    
    def modifyeq_dissolve_only(self):
        return {}
    
    def kill_phrqc_workers(self):
        self.phrqc.kill()
        return
    
class phase(object):
    """
    phase_type,phase_id, mvol,phaseqty
    """
    def __init__(self,**kwargs):
        for k,v in kwargs.items():
            setattr(self,k,v)
    
    @property
    def c(self):
        if self.type == 'diffusive':
            return self.inst._diffusive_phaseqty[self.phaseid]
        elif self.type == 'non_diffusive':
            return self.inst._non_diffusive_phaseqty[self.phaseid]

    @c.setter
    def c(self,v):
        if self.type=='diffusive':
            v = v *np.ones(self.shape)
            self.inst._diffusive_phaseqty[self.phaseid]=np.asfortranarray(v)
        elif self.type == 'non_diffusive':
            self.inst._non_diffusive_phaseqty[self.phaseid]=np.asfortranarray(v)
    
    @property
    def vol(self):
        return self.c * self.mvol
    
    @property
    def phrqc_dphase(self):
        if self.type == 'diffusive':
            return self.inst._diffusive_phrqc_dphases[self.phaseid]
        elif self.type == 'non_diffusive':
            return self.inst._non_diffusive_phrqc_dphases[self.phaseid]
    
    @phrqc_dphase.setter
    def phrqc_dphase(self,v):
        if self.type=='diffusive':
            v = v *np.ones(self.shape)
            self.inst._diffusive_phrqc_dphases[self.phaseid]=np.asfortranarray(v)
        elif self.type == 'non_diffusive':
            self.inst._non_diffusive_phrqc_dphases[self.phaseid]=np.asfortranarray(v)
    
class Solid(object):
    """
    """
    def __init__(self,domain,domain_params,solver_params):
        self.d = domain.d
        if self.d == 2:
            self.shape = (domain.ny,domain.nx)
        else:
            self.shape = (domain.nz,domain.ny,domain.nx)
        self.nodetype = domain.nodetype
        self.solid_phases = domain_params.get('solid_phases',{})
        self.voxel_vol = domain_params.get('voxel_vol', 1)
        self.frac = domain_params.get('frac',0.5)
        self._set_fmethods()
        self.nsolids = np.sum(self.nodetype>0)
        self.reassign_nodetype()

    def _set_fmethods(self):
        if self.d == 2:
            update = yantra._solvers.update2d
        elif self.d == 3:
            update = yantra._solvers.update3d
        self.reassign_nodetype= Fmethod(self,update.reassign_nodetype,['_nodetype'],)
        self.update_nodetype=Fmethod(self,update.update_nodetype,
                                     ['_nodetype','non_diff_solid_vol','voxel_vol','frac'])
        self.update_non_diffusive_phases=Fmethod(self,update.update_non_diffusive_phases,
            ['_non_diffusive_phrqc_dphases','_non_diffusive_phaseqty','_nodetype','non_diff_solid_vol','non_diff_mvol','voxel_vol'])
        self.update_diffusive_phases=Fmethod(self,update.update_diffusive_phases,
            ['_diffusive_phrqc_dphases','_diffusive_phaseqty','_nodetype','diff_solid_vol','diff_mvol'])
        self._phaseqty_for_phrqc=Fmethod(self,update.phaseqty_for_phrqc,
                                         ['nodetype','phaseqty'],inplace_update=False)
        self._poros=Fmethod(self,update.porosity,['diff_solid_vol','voxel_vol','nodetype'],inplace_update=False)                
        
    def update(self,dphases):
        """
        updates solid phases
        """
        if self.nphases >0:
            self.phrqc_dphases = dphases
            self.prev_nsolids=self.nsolids
            if self.n_non_diffusive_phases > 0:
                self.update_non_diffusive_phases()
                self.update_nodetype()
            if self.n_diffusive_phases > 0:
                self.update_diffusive_phases()
            self.nsolids=np.sum(self.nodetype>0)
            if self.nsolids != self.prev_nsolids:
                self.reassign_nodetype()
                return self.phaseqty_for_phrqc()
            else:
                return []
        else: 
            return []
    @property
    def poros(self):
        if self.n_diffusive_phases > 0:
            return self._poros()
        else:
            return 1
    def phaseqty_for_phrqc(self):
        out=self._phaseqty_for_phrqc()
        arranged_out = {}
        for i,p in enumerate(self.phaselist):
            arranged_out[p]=out[i]
        return arranged_out
    
    @property
    def nodetype(self):
        return self._nodetype
    
    @nodetype.setter
    def nodetype(self,val):
        self._nodetype = np.asfortranarray(val)
    
    @property
    def solid_phases(self):
        s_phases ={}
        for p in self.phaselist:
            phase = getattr(self,p)
            s_phases[p]={}
            s_phases[p]['type'] = phase.type
            s_phases[p]['c'] = phase.c
            s_phases[p]['mvol'] = phase.mvol
            s_phases[p]['phrqc_dphase'] = phase.phrqc_dphase
        return s_phases
        
    @solid_phases.setter
    def solid_phases(self,val):
        self.non_diffusive_phase_list=[]
        self.diffusive_phase_list=[]
        itr_diff=0
        itr_non_diff = 0
        for k,v in val.items():
            if v['type'] == 'diffusive':
                phaseid=itr_diff
                itr_diff +=1
                self.diffusive_phase_list.append(k)
            if v['type'] == 'non_diffusive':
                phaseid=itr_non_diff
                itr_non_diff +=1                
                self.non_diffusive_phase_list.append(k)
            setattr(self,k,phase(inst=self,type=v['type'],
                     mvol=v['mvol'],phaseid=phaseid,shape=self.shape))
        self.phaselist= self.diffusive_phase_list+self.non_diffusive_phase_list
        self.diffusive_phaseqty_shape = tuple([self.n_diffusive_phases]+list(self.shape))
        self.non_diffusive_phaseqty_shape = tuple([self.n_non_diffusive_phases]+list(self.shape))
        self.phaseqty_shape = tuple([self.nphases]+list(self.shape))
        self.phaseqty = val
        
    @property
    def nphases(self):
        return len(self.phaselist)
    
    @property
    def n_non_diffusive_phases(self):
        return len(self.non_diffusive_phase_list)

    @property
    def n_diffusive_phases(self):
        return len(self.diffusive_phase_list)

    @property
    def diff_mvol(self):
        try: 
            getattr(self,'_diff_mvol')
        except AttributeError:
            mvol = []
            for p in self.diffusive_phase_list:
                mvol.append(getattr(self,p).mvol)
            self._diff_mvol = np.asfortranarray(np.array(mvol))
        return self._diff_mvol

    @property
    def non_diff_mvol(self):
        try: 
            getattr(self,'_non_diff_mvol')
        except AttributeError:
            mvol = []
            for p in self.non_diffusive_phase_list:
                mvol.append(getattr(self,p).mvol)
            self._non_diff_mvol = np.asfortranarray(np.array(mvol))
        return self._non_diff_mvol
        
    @property
    def tot_solid_vol(self):
        return self.diff_solid_vol + self.non_diff_solid_vol
        
    @property
    def diff_solid_vol(self):
        try:
            getattr(self,'_diff_solid_vol')
        except AttributeError:
            self._diff_solid_vol=0
            for name in self.diffusive_phase_list:
                val=getattr(self,name)
                self._diff_solid_vol += val.mvol * val.c
        return np.asfortranarray(self._diff_solid_vol)

    @diff_solid_vol.setter
    def diff_solid_vol(self,v):
        self._diff_solid_vol = np.asfortranarray(v*np.ones(self.shape))

    @property
    def non_diff_solid_vol(self):
        try:
            getattr(self,'_non_diff_solid_vol')
        except AttributeError:
            self._non_diff_solid_vol=0
            for name in self.non_diffusive_phase_list:
                val=getattr(self,name)
                self._non_diff_solid_vol += val.mvol * val.c
        return np.asfortranarray(self._non_diff_solid_vol)
    
    @non_diff_solid_vol.setter
    def non_diff_solid_vol(self,v):
        self._non_diff_solid_vol = np.asfortranarray(v*np.ones(self.shape))


    @property        
    def phaseqty(self):
        if (self.n_non_diffusive_phases > 0) and  (self.n_diffusive_phases > 0):
            return np.concatenate((self._diffusive_phaseqty,self._non_diffusive_phaseqty),axis=0)
        elif  (self.n_non_diffusive_phases > 0) and  (self.n_diffusive_phases == 0):
            return self._non_diffusive_phaseqty
        elif  (self.n_non_diffusive_phases == 0) and  (self.n_diffusive_phases > 0):
            return self._diffusive_phaseqty
        
        
    @phaseqty.setter
    def phaseqty(self,val):
        """
        sets phase quantity
        """
        try:
            getattr(self,'_diffusive_phaseqty')
        except AttributeError:
            self._diffusive_phaseqty = np.asfortranarray(np.zeros(self.diffusive_phaseqty_shape))
        try:
            getattr(self,'_non_diffusive_phaseqty')
        except AttributeError:
            self._non_diffusive_phaseqty = np.asfortranarray(np.zeros(self.non_diffusive_phaseqty_shape))
        if type(val) is dict:
            for k,v in val.items():
                    phase=getattr(self,k)
                    if type(v) is dict:
                        phase.c = v['c']
                    else:
                        phase.c = v
        else:
            assert val.shape(0)==self.nphases
            for i,p in enumerate(self.phaselist):
                phase=getattr(self,p)
                phase.c=v[i]
            
    @property
    def phrqc_dphases(self):
        if (self.n_non_diffusive_phases > 0) and  (self.n_diffusive_phases > 0):
            return np.concatenate((self._diffusive_phrqc_dphases,self._non_diffusive_phrqc_dphases),axis=0)
        elif  (self.n_non_diffusive_phases > 0) and  (self.n_diffusive_phases == 0):
            return self._non_diffusive_phrqc_dphases
        elif  (self.n_non_diffusive_phases == 0) and  (self.n_diffusive_phases > 0):
            return self._diffusive_phrqc_dphases
    
    @phrqc_dphases.setter
    def phrqc_dphases(self,val):
        try:
            getattr(self,'_non_diffusive_phrqc_dphases')
        except AttributeError:
            self._non_diffusive_phrqc_dphases = np.asfortranarray(np.zeros(self.non_diffusive_phaseqty_shape))
        try:
            getattr(self,'_diffusive_phrqc_dphases')
        except AttributeError:
            self._diffusive_phrqc_dphases = np.asfortranarray(np.zeros(self.diffusive_phaseqty_shape))
        if type(val) is dict:
            for k,v in val.items():
                    phase=getattr(self,k)
                    if type(v) is dict:
                        phase.phrqc_dphase = v['c']
                    else:
                        phase.phrqc_dphase = v
        else:
            assert val.shape[0]==self.nphases
            for i,p in enumerate(self.phaselist):
                phase=getattr(self,p)
                phase.phrqc_dphase=val[i]
            

        