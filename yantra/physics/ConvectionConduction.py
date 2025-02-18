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
from yantra import __version__,__license__
__name__ = 'ConvectionConduction'
__author__ = 'Ravi Patel'   
__email__ = 'ravee.a.patel@gmail.com'
        
class ConvectionConduction(LB, metaclass=LBMeta):
    """
    Solves convection conduction equation written as
    ..math::
        \partial_{t}T + u  \cdot \nabla T = \frac{1}{\rho c_p} \nabla \kappa \cdot \nabla T + ss
    Flux
    ..math::
        \vec{j}=\rho c_p u T -  \kappa \cdot \nabla T
    where,
    T = temperature ..math::[\theta^{1}]
    u = velocity ..math::[L^{1}T^{-1}] 
    ..math::c_p = specific heat capacity at constant pressure ..math::[L^2T^{-2}\theta^{-1}]
    ..math::\rho = density of the material(fluid/solid) ..math::[M^1L^-3]
    ..math:: \kappa = thermal conductivity ..math::[M^1L^1T^{-3}\theta^{-1}]
    j = total heat flux ..math::[N^{1}T^{-3}]
    ss = Source/sink term ..math::[\theta^{1}T^{-1}]
    The bounce back condition in streaming step imposes following condition for all nodetype > 0
    {\nabla}T  \cdot {n} = 0 
    """
    _signature = 'yantra.physics.AdvectionDiffusion.AdvectionDiffusion'
    #list of variables for ccd class 
    _vars = {'T':Variable('scalar',0,{'theta':1},True,False,'domain_params','concentration'),
             'u':Variable('vector',0,{'L':1,'T0':-1},True,False,'domain_params','velocity'),
             'rho':Variable('scalar',1.,{'M':1,'L':-3},True,False,'domain_params','density of the phases'),
             'cp':Variable('scalar',1.,{'L':2,'T0':-2,'theta':-1},True,False,'domain_params','specific heat capacity'),
             'kappa':Variable('scalar',1./6.,{'M':1,'L':1,'T0':-3,'theta':-1},True,False,'domain_params','diffusion coefficient'),
             'alpha_ref':Variable('parameter',1./6.,{'M':1,'L':1,'T0':-3,'theta':-1},True,False,'solver_params','reference diffusion coefficient for setting timestep'),
             'ss':Variable('scalar',0,{'theta':1,'T0':-1},True,False,'domain_params','source/sink term'),
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
             'tfactbased': Variable('flag',0.,{},False,False,'solver_params','to set time scale according to factor of tfact*dx**2/alpha'),
             'tfact': Variable('parameter',1./6.,{},False,False,'solver_params','to set time scale according to factor of tfact*dx**2/alpha'),
             'lattice': Variable('parameter','D2Q5',{},False,False,'solver_params','lattice type for the model'),
             'f': Variable('dist_func',0,{},False,True,'domain_params','LB distribution function'),
             'time': Variable('parameter',0,{},False,True,'solver_params','simulation time in phyiscal units'),
             'iters': Variable('parameter',0,{},False,True,'solver_params','number of iterations in simulations'),
            }
            
    def __init__(self,domain,domain_params,bc_params,solver_params):
        """
        Initializes AdvectionDiffusion equation selfance
        """
        domain=self._convert_domaintype(domain)
        self._construct_fort_solver(solver_params)
        #distribution function default value as equilibrium distribution function
        self._vars['f'].default=self.compute_edf
        #call LB initialization method
        LB.__init__(self,domain,domain_params,bc_params,solver_params)
        self._set_relaxation_params()
        
    def _set_relaxation_params(self):
        """
        set relaxation parameter
        """
        self.tau = self.diff2tau(alpha=self._alpha)
            
    def _get_conv_factors(self,*args):
        """
        computes conversion factors for different variables
        """
        if len(args)==3:
            domain,domain_params,solver_params=args
        klist=['rho','cp','kappa','tauref','tfact','tfactbased']
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
        rho,cp,kappa,tauref,tfact,tfactbased = vals
        alpha = kappa/(rho*cp)
        alpha_ref = np.max(alpha)
        if len(args)==3: alpha_ref = solver_params.get('alpha_ref',alpha_ref)
        if len(args)==3:
            self._vars['alpha_ref'].default = alpha_ref
        else:
            self.alpha_ref = alpha_ref
        convfactors={}
        if tfactbased:
            dt= tfact * dx**2/alpha_ref
        else:
            alpha_ref_lb=self.tau2diff(tau = tauref)
            alpha0= alpha_ref/alpha_ref_lb
            dt=dx**2/alpha0
        base = {}
        base['L']= dx
        base['M']=1
        base['theta']=1        
        base['T0']= dt
        convfactors = deepcopy(base)
        for k,v in self._vars.items():
            if v.isphyvar:
                cfact = 1
                for p,v in v.dimension.items():
                    cfact*=base.get(p,1)**v  
                convfactors[k]=cfact
        self.convfactors= convfactors
        self.convfactors['tot_flux']=base['M']*base['T0']**(-3)
        self.convfactors['grad_T']=base['theta']*base['L']**(-1)
                
    def tau2diff(self,**kwargs):
        """
        Gets  diffusion coefficient(in LB units) from tau
        
        Parameters
        ----------
        tau: int, float or ndarray
            relaxation parameter
        
        Returns
        -------
        alpha:    int, float or ndarray
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
        alpha: int,float or ndarray
            in LB units

        Returns
        -------
        tau:    int, float or ndarray
            relaxation parameter             
        """
        if 'alpha' in kwargs:
            alpha = kwargs['alpha']
        else:
            alpha = self._alpha
        es2=getattr(self,'es2',self._vars['es2'].default)
        return alpha/es2+0.5

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
            ccd = getattr(_solvers,'ccd2d')
        elif d==3:
            ccd = getattr(_solvers,'ccd3d')   
        #compute_macro_var
        self.compute_macro_var = Fmethod(self,ccd.compute_macro_var,
                                        ['_f', '_T','nodetype'])
        #compute edf
        self.compute_edf = Fmethod(self,ccd.compute_edf,
                                        ['_T', '_u','nodetype'],inplace_update=False)
        #collison
        if collision_model == 'srt':
            self.collide = Fmethod(self,ccd.collide_srt,
                                   ['_f', '_T',  '_u','nodetype','_tau','_ss','_kappa', '_rho','_cp','_periodicity'])
            
        if collision_model == 'trt':
            self.collide = Fmethod(self,ccd.collide_trt,
                                   ['_f', '_T',  '_u','nodetype','_tau','_magic_para','_ss','_kappa', '_rho','_cp','_periodicity'])      
        #stream
        self.stream = Fmethod(self,ccd.stream_and_bounce,['_f', 'nodetype','_periodicity'])      
        #apply_bc
        self.apply_bc = Fmethod(self,ccd.boundary_conditions,
        ['_f', '_u', 'nodetype', '_tau','_kappa', '_rho','_cp', '_interp'],['_bc'])
        #total_heat_flux
        self._get_tot_flux = Fmethod(self,ccd.total_flux,
        ['_f', '_u', '_rho','_cp','_kappa', '_tau', 'nodetype'])
        #gradient
        self._get_grad_t = Fmethod(self,ccd.grad_t,
        ['_f', '_u', '_tau', 'nodetype'])        
        
        
    def advance(self):
        self.time +=self.dt
        self.iters+=1
        self.collide()
        self.stream()
        self.apply_bc()
        self.compute_macro_var()
    
    def steady_state(self,verbose=True,check_iters=100,verbose_iters=500,tol =1e-6,
                     err_cal_method='l2norm',skip_boundaries=True):
        while 1:
            T_old = self._T.copy()
            if skip_boundaries:
                if self.d==2:
                    T_old = T_old[1:-1,1:-1]
                elif self.d==3:
                    T_old = T_old[1:-1,1:-1,1:-1]  
            self.run(iters=self.iters+check_iters)
            T = self._T
            if skip_boundaries:
                if self.d==2:
                    T = T[1:-1,1:-1]
                elif self.d==3:
                    T = T[1:-1,1:-1,1:-1] 
            if self.iters>1:
                err=getattr(self,err_cal_method)(T, T_old)
                if (err/tol) <1: 
                    break
                if verbose and (self.iters%verbose_iters==0):
                    print('Iteration No:%s, Error: %s'%(self.iters,err))    
    @property
    def alpha(self):
        return self.kappa/(self.rho*self.cp)
    
    @property
    def _alpha(self):
        return self._kappa/(self._rho*self._cp)
    
    @property
    def _grad_T(self):
        return self._get_grad_t()
    
    @property
    def grad_T(self):
        return self._grad_T * self.convfactors['grad_T']

    @property
    def _tot_flux(self):
        return self._get_tot_flux()
    
    @property
    def tot_flux(self):
        return self._tot_flux * self.convfactors['tot_flux']
    
    @property
    def dt(self):
        """
        """
        return self.convfactors['T0']
    
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
            if name=='T':
                val=self._bc[i+'val']*self.convfactors['T']       
            elif self._bc[i+'bc']=='tot_flux':
                val=self._bc[i+'val']*self.convfactors['tot_flux'] 
            elif self._bc[i+'bc']=='grad_T':
                val=self._bc[i+'val']*self.convfactors['grad_T'] 
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
            if self.d == 2:
                if i in ['left','right']:
                    ones=np.array([1]*self.ny)
                elif i in ['top','bottom']:
                    ones=np.array([1]*self.nx)
            elif self.d==3:
                if i in ['left','right']:
                    ones=np.ones(self.nz,self.ny)
                elif i in ['top','bottom']:
                    ones=np.ones(self.nz,self.nx)
                elif i in ['front','back']:
                    ones=np.ones(self.ny,self.nx)
            if self._bc[i+'bc']=='T':
                self._bc[i+'val']=bc[i][1]/self.convfactors['T']*ones       
            elif self._bc[i+'bc']=='tot_flux':
                self._bc[i+'val']=bc[i][1]/self.convfactors['tot_flux'] * ones
            elif self._bc[i+'bc']=='grad_T':
                self._bc[i+'val']=bc[i][1]/self.convfactors['grad_T'] * ones        
            else:
                self._bc[i+'val']=bc[i][1] * ones
                

