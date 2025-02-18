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

from yantra import _solvers
from yantra._base import LBMeta,Fmethod,Variable
from copy import deepcopy
from .AdvectionDiffusion import AdvectionDiffusion
import numpy as np
from copy import deepcopy as copy
from yantra import __version__,__license__
__name__ = 'Diffusion'
__author__ = 'Ravi Patel'   
__email__ = 'ravee.a.patel@gmail.com'

class MultilevelAdvectionDiffusion(AdvectionDiffusion, metaclass=LBMeta):
    """
    Solves advection diffusion equation written as
    ..math::
        \partial_{t} \phi c = -\vec{nabla}  \cdot \vec{j}+ ss
        \vec{j}=  \phi c\vec{u} - \phi \zeta_{a} D_{0}\vec{\nabla}c 
    where,
    c = concentration in the aqeuous phase ..math::[N^{1}L^{-3}]
    u = velocity of the fluid (not darcy velocity) ..math::[L^{1}T^{-2}]
    j = flux ..math::[N^{1}L^{-2}T^{-1}]
    $\phi$ = porosity of the media
    $\zeta_{a}$ = apparent tortousity of the media [-]
    $D_{0}$ = diffusion coefficient in pore water..math::[L^{2}T^{-1}]
    ss = Source/sink term ..math::[N^{1}L^{-3}T^{-1}]
    The term \phi \zeta_{a} D_{0} is referred to as effective diffusivity (D_e)
    
    The bounce back condition in streaming step imposes following condition for all nodetype > 0
    \vec{\nabla}c  \cdot \hat{n} = 0 
    """
    _signature = 'yantra.physics.MultilevelAdvectionDiffusion.MultilevelAdvectionDiffusion'
    _vars=deepcopy(AdvectionDiffusion._vars)
    _vars['D0']= _vars.pop('D')
    _vars['D0'].type = 'parameter'
    _vars['Deref']=_vars.pop('Dref')
    _vars['collision_model'].default='trt'
    _vars['poros']=Variable('scalar',1.,{},False,False,'domain_params','porosity')
    _vars['dporosdt']=Variable('scalar',0.,{'T':-1},True,False,'domain_params','derivative of porosity with time')
    _vars['dcdt']=Variable('scalar',0.,{'N':1,'L':-3,'T':-1},True,False,'domain_params','derivative of porosity with time')
    _vars['cphi_fact']=Variable('parameter',1./3.5,{},False,False,'solver_params','cphi fact for TRT multilevel formulation')
    _vars['cphi']=Variable('parameter',1./3.5,{},False,False,'solver_params','cphi_fact*max(poros)')
    _vars['const_cphi']=Variable('parameter',0,{},False,False,'solver_params','cphi_fact*max(poros)')
    _vars['app_tort']= Variable('scalar',1.,{},False,False,'domain_params','apparent tortousity')
    #params related to diffuse velocity    
#    _vars.pop('u_adv')
#    _vars.pop('Dref')
        
    def _set_relaxation_params(self):
        """
        sets relaxation parameter based on effective diffusion
        """
        collision_model = self.collision_model
        if collision_model.lower()=='diff_vel':
            self.tau = self.diff2tau(De=self._Deref)
        else:
            self.tau = self.diff2tau(De=self._De)       

    def _get_conv_factors(self,*args):
        """
        computes conversion factors for different variables
        """
#TODO: correct unit conversion for trt model
        if len(args)==3:
            domain,domain_params,solver_params = args
        klist=['D0','poros','app_tort','cphi_fact','cphi','const_cphi','tauref','tfactbased','tfact','collision_model','d']
        vals=[] 
        if len(args)==3:
            for k in klist:
                if self._vars[k].belongs_to=='domain_params':
                   vals.append(domain_params.get(k,self._vars[k].default)) 
                elif self._vars[k].belongs_to=='solver_params':
                   vals.append(solver_params.get(k,self._vars[k].default))    
            dx = getattr(domain,'dx')
        else:
            for k in klist:
                vals.append(getattr(self,k))
            dx = self.dx
            De = self.De
        D0,poros,app_tort,cphi_fact,cphi,const_cphi,tauref,tfactbased,tfact,collision_model,d=vals
        if len(args)==3: De=(D0*poros*app_tort*np.ones(domain.shape))
        De[De==0]=np.NaN
        Deref= np.nanmax(De)
        if len(args)==3: Deref= solver_params.get('Deref',Deref) 
        poros = deepcopy(poros)
        if len(args)==3: poros = poros*np.ones(domain.shape)
        poros[poros==0]=np.NaN
        collision_model = collision_model.lower()
        if (collision_model == 'trt') and ( not const_cphi):
            cphi =  np.nanmin(poros)*cphi_fact
        elif collision_model == 'diff_vel':
            if d == 2 : cphi = 1./3.
            if d == 3: cphi = 1./3.5
        if len(args)==3:
            self._vars['cphi'].default = cphi
            self._vars['Deref'].default= Deref
        else:
            self.cphi = cphi
            self.Deref= Deref
        convfactors={}
        if tfactbased:
            dt= tfact * dx**2/Deref
        else:
            Dereflb=self.tau2diff(tau = tauref,cphi=cphi)
            D0= Deref/Dereflb
            dt=dx**2/D0
        base = {}
        base['L']= dx
        base['N']=base['L']**3
        base['T']= dt
        convfactors = deepcopy(base)
        for k,v in self._vars.items():
            if v.isphyvar:
                cfact = 1
                for p,v in v.dimension.items():
                    cfact*=base.get(p,1)**v  
                convfactors[k]=cfact
        self.convfactors= convfactors
    
    def _update_transport_params_trt(self,poros,app_tort,auto_time_step=True):
        self.poros = poros
        self.app_tort = app_tort
        if auto_time_step:
            bc = deepcopy(self.bc)
            old_convfactors = deepcopy(self.convfactors)
            self._get_conv_factors()
            #reset values in LB units
            for key,val in self._vars.items():
                if val.isphyvar:
                    try:
                        v=getattr(self,'_'+key)
                        v*=(old_convfactors[key]/self.convfactors[key])
                        setattr(self,'_'+key,v)          
                    except AttributeError:
                        pass
            self.bc = bc
        self._set_relaxation_params()

    
    def _update_transport_params_diff_vel(self,poros,app_tort,auto_time_step=True):
        self.poros = poros
        self.app_tort = app_tort
        self.Dr = self.De - self.Deref
#        if auto_time_step:
#            bc = deepcopy(self.bc)
#            old_convfactors = deepcopy(self.convfactors)
#            self._get_conv_factors()
#            #reset values in LB units
#            for key,val in self._vars.iteritems():
#                if val.isphyvar:
#                    try:
#                        v=getattr(self,'_'+key)
#                        v*=(old_convfactors[key]/self.convfactors[key])
#                        setattr(self,'_'+key,v)          
#                    except AttributeError:
#                        pass
#            self.bc = bc
#        self._set_relaxation_params()

    def update_transport_params(self,poros,app_tort,auto_time_step=True):
        collision_model = self.collision_model.lower()
        if collision_model == 'diff_vel':
           self._update_transport_params_diff_vel(poros,app_tort,auto_time_step=True)
        elif collision_model == 'trt':
            self._update_transport_params_trt(poros,app_tort,auto_time_step=True)
        return
        
    @property
    def De(self):
        """
        effective diffusivity in physical units
        """
        return self.D0 * self.app_tort * self.poros
    
    @property
    def _De(self):
        """
        effective diffusivity in physical units
        """
        return self._D0 * self.app_tort * self.poros


    def tau2diff(self,**kwargs):
        """
        Gets  diffusion coefficient(in LB units) from tau
        
        Parameters
        ----------
        tau: int, float or ndarray
            relaxation parameter
        
        Returns
        -------
        D:    int, float or ndarray
            Diffusion coefficient in LB units
        """
        key_args=['tau','cphi']
        val =[]
        for k in key_args:
            try:
                val.append(kwargs[k])
            except KeyError:
                val.append(kwargs.get(k,getattr(self,k)))           
        return val[1]*(val[0]-0.5)
        
    def diff2tau(self,**kwargs):
        """
        gets tau from diffusion coefficient (in LB units)

        Parameters
        ----------
        De: int,float or ndarray
            effective diffusion in LB units

        Returns
        -------
        tau:    int, float or ndarray
            relaxation parameter             
        """
        key_args=['De','cphi']
        val =[]
        for k in key_args:
            try:
                val.append(kwargs[k])
            except KeyError:
                if k=='De': k ='_'+k
                val.append(kwargs.get(k,getattr(self,k)))
        return val[0]/val[1]+0.5

    
    def _construct_fort_solver(self,solver_params):
        """
        constructs fortran solver based on given inputs
        """
        d = solver_params.get('d',self._vars['d'].default)
        lattice = solver_params.get('lattice',self._vars['lattice'].default)
        collision_model = solver_params.get('collision_model',self._vars['collision_model'].default)
        lattice = lattice.upper()
        collision_model = collision_model.lower()
        if int(lattice[1]) == 2 and d == 2:
                try:
                    assert lattice == 'D2Q5'
                except AssertionError:
                    raise ValueError("Only D2Q5 lattice implementation exists for 2D")
        elif int(lattice[1]) == 3 and d == 2:
            self._vars['d'].default=d=3
            try:
                assert lattice == 'D3Q7'
                self._vars['q'].default=7
            except AssertionError:
                raise ValueError("Only D3Q7 lattice implementation exists for 3D")
        elif int(lattice[1]) == 3 and d == 3:
            try:
                assert lattice == 'D3Q7'
                self._vars['q'].default=7
            except AssertionError:
                raise ValueError("Only D3Q7 lattice implementation exists for 3D")
        elif int(lattice[1]) == 2 and d == 3:
            self._vars['lattice'].default=lattice='D3Q7'
            self._vars['q'].default=7
        if collision_model == 'trt':
            if d==2:
                ade = getattr(_solvers,'multilevel_ade2d')
            elif d==3:
                ade = getattr(_solvers,'multilevel_ade3d')   
            #compute_macro_var
            self.compute_macro_var = Fmethod(self,ade.compute_macro_var,
                                            ['_f', '_c', '_flux', '_u','_poros','nodetype','_tau'])
            #compute edf
            self.compute_edf = Fmethod(self,ade.compute_edf,
                                            ['_c', '_u','_cphi','_poros','nodetype'],inplace_update=False)
    
            #collison
            self.collide = Fmethod(self,ade.collide,
                               ['_f', '_c',  '_u','_cphi','_poros','nodetype',
                                '_tau','_magic_para','_ss'])      
            #stream
            self.stream = Fmethod(self,ade.stream_and_bounce,['_f', 'nodetype','_periodicity'])      
            #apply_bc
            self.apply_bc = Fmethod(self,ade.boundary_conditions,
            ['_f', '_u', 'nodetype', '_poros','_cphi','_tau', '_interp','grid_type'],['_bc'])
        elif collision_model == 'diff_vel':
            if d==2:
                ade = getattr(_solvers,'ade2d')
            elif d==3:
                ade = getattr(_solvers,'ade3d')   
            #compute_macro_var
            self.compute_macro_var = Fmethod(self,ade.compute_macro_var,
                                            ['_f', '_c', '_flux', '_u','nodetype','_tau'])
            #compute edf
            self.compute_edf = Fmethod(self,ade.compute_edf,
                                            ['_c', '_u','nodetype'],inplace_update=False)
            #collison
            self.collide = Fmethod(self,ade.collide_diff_vel,
                                   ['_f', '_c', '_Dr', '_u','_u_adv','nodetype','_tau','_ss'])
            #stream
            self.stream = Fmethod(self,ade.stream_and_bounce,['_f', 'nodetype','_periodicity'])      
            #apply_bc
            self.apply_bc = Fmethod(self,ade.boundary_conditions,
            ['_f', '_u', 'nodetype', '_tau', '_interp','grid_type'],['_bc'])           
        else:
            raise ValueError("%s collision model not available"%collision_model)
 
    
    def _compute_ss_poros(self):
        self._ss+= (1.-self._poros)*self._dcdt  - self._c * self._dporosdt

    def _advance_diff_vel(self):
        c_old = copy(self.c)
        self.compute_macro_var()
        dcdt =(self.c-c_old)/self.dt
        self.dcdt = (self.dcdt+dcdt)/2.
        self._compute_ss_poros()
        self.collide()
        self.stream()
        self.apply_bc()
        self.ss = 0
    
    def _advance_trt(self):
        self.compute_macro_var()
        self.collide()
        self.stream()
        self.apply_bc()
    
    def advance(self):
        collision_model = self.collision_model.lower()
        if collision_model == 'diff_vel':
           self._advance_diff_vel()
        elif collision_model == 'trt':
            self._advance_trt()
        self.time +=self.dt
        self.iters +=1
        return
    
        