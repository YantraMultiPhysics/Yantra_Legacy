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
#
#=======================================================================================

from copy import deepcopy  
from collections import defaultdict
from yantra import _solvers
import numpy as np
from yantra._base import LB, LBMeta,Variable,Fmethod
import warnings
from yantra import __version__,__license__
__name__ = 'AdvectionDiffusion'
__author__ = 'Ravi Patel'   
__email__ = 'ravee.a.patel@gmail.com'
        
class AdvectionDiffusion(LB, metaclass=LBMeta):
    """
    Solves advection diffusion equation written as
    ..math::
        \partial_{t}c = -\vec{nabla} \cdot \vec{j}+ ss
        \vec{j}=  c\vec{u} - D \vec{\nabla}c 
    where,
    c = concentration ..math::[N^{1}L^{-3}]
    u = velocity ..math::[L^{1}T^{-2}]
    j = flux ..math::[N^{1}L^{-2}T^{-1}]
    D = diffusion coefficient ..math::[L^{2}T^{-1}]
    ss = Source/sink term ..math::[N^{1}L^{-3}T^{-1}]
    The bounce back condition in streaming step imposes following condition for all nodetype > 0
    \vec{\nabla}c  \cdot \hat{n} = 0 
    """
    _signature = 'yantra.physics.AdvectionDiffusion.AdvectionDiffusion'
    #list of variables for ADE class 
    _vars = {'c':Variable('scalar',0,{'N':1,'L':-3},True,False,'domain_params','concentration'),
             'flux':Variable('vector',0,{'N':1,'L':-2,'T':-1},True,False,'domain_params','flux'),
             'u':Variable('vector',0,{'L':1,'T':-1},True,False,'domain_params','velocity'),
             'D':Variable('scalar',1./6.,{'L':2,'T':-1},True,False,'domain_params','diffusion coefficient'),
             'Dref':Variable('parameter',1./6.,{'L':2,'T':-1},True,False,'solver_params','reference diffusion coefficient for setting timestep'),
             'ss':Variable('scalar',0,{'N':1,'L':-3,'T':-1},True,False,'domain_params','source/sink term'),
             'tau':Variable('scalar',1,{},False,False,'solver_params',
                            'relaxation parameter/for TRT scheme its anti-symmetric component'),
             'tauref':Variable('parameter',1,{},False,False,'solver_params',
                            'refrence relaxation parameter/for TRT scheme its anti-symmetric component which is used to set timestep'),
             'q':Variable('int',5,{},False,False,'solver_params','Number of lattice directions'),
             'd': Variable('int',2,{},False,False,'solver_params','Dimension of the domain'),
             'es2':Variable('parameter',1./3.,{},False,False,'solver_params','pseudo velocity of the sound'),
             'magic_para':Variable('parameter',1./6.,{},False,False,'solver_params','Magic parameter for TRT LBM model'),
             'collision_model': Variable('parameter','SRT',{},False,False,'solver_params','model for LB collision term'),
             'interp': Variable('flag',1,{},False,False,'solver_params','interpolate to set correct boundary'),
             'tfactbased': Variable('flag',0.,{},False,False,'solver_params','to set time scale according to factor of tfact*dx**2/D'),
             'tfact': Variable('parameter',1./6.,{},False,False,'solver_params','to set time scale according to factor of tfact*dx**2/D'),
             'lattice': Variable('parameter','D2Q5',{},False,False,'solver_params','lattice type for the model'),
             'Dr': Variable('scalar',0.,{'L':2,'T':-1},True,False,'domain_params','remaing part of diffusion coefficient for diffusion velocity formulation'),
             'u_adv': Variable('vector',0.,{'L':1,'T':-1},True,False,'domain_params','advective velocity specified in ade'),
             'f': Variable('dist_func',0,{},False,True,'domain_params','LB distribution function'),
             'time': Variable('parameter',0,{},False,True,'solver_params','simulation time in phyiscal units'),
             'iters': Variable('parameter',0,{},False,True,'solver_params','number of iterations in simulations'),
            }
            
    def __init__(self,domain,domain_params,bc_params,solver_params):
        """
        Initializes AdvectionDiffusion equation instance
        """
        domain=self._convert_domaintype(domain)
        self._construct_fort_solver(solver_params)
        #distribution function default value as equilibrium distribution function
        self._vars['f'].default=self.compute_edf
        #call LB initialization method
        LB.__init__(self,domain,domain_params,bc_params,solver_params)
        #adjustment to be made for diffusion velocity formulation
        collision_model = self.collision_model
        if collision_model.lower()=='diff_vel':
            try:
                self.Dr = self.D - self.Dref
                self.u_adv = deepcopy(self.u)
            except KeyError:
                #Added to make it compatible with multilevel class
                self.Dr = self.De - self.Deref
                self.u_adv = deepcopy(self.u)                
#        else:
#            try:
#                del self.Dr,self.u_adv
#            except: AttributeError
        self._set_relaxation_params()
        
    def _set_relaxation_params(self):
        """
        set relaxation parameter
        """
        collision_model = self.collision_model
        if collision_model.lower()=='diff_vel':
            self.tau = self.diff2tau(D=self._Dref)
        else:
            self.tau = self.diff2tau(D=self._D)
            
    def _get_conv_factors(self,*args):
        """
        computes conversion factors for different variables
        """
        if len(args)==3:
            domain,domain_params,solver_params=args
        klist=['D','tauref','tfact','tfactbased']
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
        D,tauref,tfact,tfactbased = vals
        Dref = np.max(D)
        if len(args)==3: Dref = solver_params.get('Dref',Dref)
        if len(args)==3:
            self._vars['Dref'].default = Dref
        else:
            self.Dref = Dref
        convfactors={}
        if tfactbased:
            dt= tfact * dx**2/Dref
        else:
            Dreflb=self.tau2diff(tau = tauref)
            D0= Dref/Dreflb
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
        if 'tau' in kwargs:
            tau = kwargs['tau']
        else:
            tau = self.tau
        es2=getattr(self,'es2',self._vars['es2'].default)
        return es2*(tau-0.5)
        
    def diff2tau(self,**kwargs):
        """
        gets tau from diffusion coefficient (in LB units)

        Parameters
        ----------
        D: int,float or ndarray
            in LB units

        Returns
        -------
        tau:    int, float or ndarray
            relaxation parameter             
        """
        if 'D' in kwargs:
            D = kwargs['D']
        else:
            D = self._D
        es2=getattr(self,'es2',self._vars['es2'].default)
        return D/es2+0.5
        
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
        if lattice=='D3Q7':
            self._vars['es2'].default=1/3.5
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
        if collision_model == 'srt':
            self.collide = Fmethod(self,ade.collide_srt,
                                   ['_f', '_c',  '_u','nodetype','_tau','_ss'])
        if collision_model == 'diff_vel':
            self.collide = Fmethod(self,ade.collide_diff_vel,
                                   ['_f', '_c', '_Dr', '_u','_u_adv','nodetype','_tau','_ss'])

        if collision_model == 'trt':
            self.collide = Fmethod(self,ade.collide_trt,
                                   ['_f', '_c',  '_u','nodetype','_tau','_magic_para','_ss'])      
        #stream
        self.stream = Fmethod(self,ade.stream_and_bounce,['_f', 'nodetype','_periodicity'])      
        #apply_bc
        self.apply_bc = Fmethod(self,ade.boundary_conditions,
        ['_f', '_u', 'nodetype', '_tau', '_interp','grid_type'],['_bc'])
    
    def advance(self):
        self.time +=self.dt
        self.iters+=1
        self.collide()
        self.stream()
        self.apply_bc()
        self.compute_macro_var()

    def steady_state(self,verbose=True,tol =1e-6,check_iters= 100,err_cal_method='l2norm',max_iters = 5000,skip_boundaries=True):
        while 1:
            c_old = self._c.copy()
            if skip_boundaries:
                if self.d==2:
                    c_old = c_old[1:-1,1:-1]
                elif self.d==3:
                    c_old = c_old[1:-1,1:-1,1:-1]  
            self.run(iters = self.iters+check_iters)
            c  = self._c.copy()
            if skip_boundaries:
                if self.d==2:
                    c = c[1:-1,1:-1]
                elif self.d==3:
                    c = c[1:-1,1:-1,1:-1]  
            err=getattr(self,err_cal_method)(c, c_old)
            if (err/tol) <1: 
                break
            if self.iters > max_iters:
                warnings.warn("solution did not converged maximum of %s iterations reached"%max_iters)
                break
            if verbose:
                print('Iteration No:%s, Error: %s'%(self.iters,err))
                
    @property
    def dt(self):
        """
        """
        return self.convfactors['T']
    @property
    def bc(self):
        """
        boundary conditions variable
        """
        bcparams = {}
        if self.d == 2:
            dirs = ['left','right','top','bottom']
        elif self.d == 3:
            dirs = ['left','right','top','bottom','front','back']
        for i in dirs:
            name = self._bc[i+'bc']
            if name=='c':
                val=self._bc[i+'val']*self.convfactors['c']       
            elif self._bc[i+'bc']=='flux':
                val=self._bc[i+'val']*self.convfactors['flux'] 
            else:
                val = self._bc[i+'val']
            bcparams[i]=[name,val]
        return bcparams
        
    @bc.setter
    def bc(self,bc_params):
        """
        setter for bc
        """
        bc=defaultdict(lambda:['nothing',0])
        if hasattr(self,'_bc'):
            bc.update(self.bc)
        bc.update(deepcopy(bc_params))
        if self.d == 2:
            dirs = ['left','right','top','bottom']
        elif self.d == 3:
            dirs = ['left','right','top','bottom','front','back']
        self._bc={}
        for name,val in bc.items():
            if (len(val)==1) and (name in dirs):val.append(0)
        for i in dirs:
            self._bc[i+'bc']=bc[i][0]  
            if self._bc[i+'bc']=='c':
                self._bc[i+'val']=bc[i][1]/self.convfactors['c']       
            elif self._bc[i+'bc']=='flux':
                self._bc[i+'val']=bc[i][1]/self.convfactors['flux'] 
            else:
                self._bc[i+'val']=bc[i][1]
                
