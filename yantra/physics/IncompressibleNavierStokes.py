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
import warnings
from yantra._base import LB, LBMeta,Variable,Fmethod
from yantra import __version__,__license__
__name__ = 'IncompressibleNavierStokes'
__author__ = 'Ravi Patel'   
__email__ = 'ravee.a.patel@gmail.com'
        
class IncompressibleNavierStokes(LB, metaclass=LBMeta):
    """
    Solves incompressible Navier-Stokes equation given by following
    ..math::
        \vec{\nabla} \cdot \vec{u} = 0 (mass conservation)
        \partial_t \vec{u} (\vec{u}\cdot\vec{\nabla}) = -(1/\rho_0) \vec{\nabla} \cdot p\hat{I}+\neu\nabla^2\vec{u} + F_{v}(momentum conservation)

    where,
    u = velocity ..math::[L^{1}T^{-1}]
    \rho_0 = density of fluid ..math::[M^{1}L^{-3}]
    p = hydrostatic pressure ..math::[M^{1}L^{-1}T^[-2]]
    \neu = Kinematic viscosity ..math::[L^{2}T^{-1}]
    F_{v} = volumetric force ..math::[M^{1}L^{-2}T^[-2]]
    The bounce back condition in streaming step imposes following condition for all nodetype > 0
    \vec{u} = 0 
    """
    _signature = 'yantra.physics.IncompressibleNavierStokes.IncompressibleNavierStokes'
    #list of variables for NS class 
    _vars = {'rhof':Variable('parameter',1,{'M':1,'L':-3},True,False,'domain_params','fluid density'),
             'rho0':Variable('parameter',1,{},True,False,'domain_params','reference density corresponding to p0 in LB units'),
             'rho':Variable('scalar',1,{'M':1,'L':-3},True,False,'domain_params','density computed from distribution function'),
             'u':Variable('vector',0,{'L':1,'T':-1},True,False,'domain_params','velocity'),
             'visc':Variable('scalar',1./6.,{'L':2,'T':-1},True,False,'domain_params','kinematic viscosity'),
             'prel':Variable('scalar',0.,{'M':1,'L':-1,'T':-2},True,False,'domain_params','relative pressure'),
             'p0':Variable('parameter',0.,{},True,False,'domain_params','reference pressure'),
             'Fv':Variable('vector',0,{'M':1,'L':-2,'T':-2},True,False,'domain_params','volumetric body force'),
             'tau':Variable('scalar',1,{},False,False,'solver_params',
                            'relaxation parameter/for TRT scheme its symmetric component'),
             'tauref':Variable('parameter',1,{},False,False,'solver_params',
                            'refrence relaxation parameter/for TRT scheme its symmetric component which is used to set timestep'),
             'q':Variable('int',9,{},False,False,'solver_params','Number of lattice directions'),
             'd': Variable('int',2,{},False,False,'solver_params','Dimension of the domain'),
             'es2':Variable('parameter',1./3.,{},False,False,'solver_params','pseudo velocity of the sound'),
             'magic_para':Variable('parameter',3./16.,{},False,False,'solver_params','Magic parameter for TRT LBM model'),
             'collision_model': Variable('parameter','SRT',{},False,False,'solver_params','model for LB collision term'),
             'forcing_model': Variable('parameter','guo',{},False,False,'solver_params','model for LB forcing term'), 
             'interp': Variable('flag',1,{},False,False,'solver_params','interpolate to set correct boundary'),
             'VelBased': Variable('flag',0.,{},False,False,'solver_params','to set time scale according to ratio of LB velocity'),
             'VelRatio': Variable('parameter',0.01,{},False,False,'solver_params','to set time scale according to factor of tfact*dx**2/D'),
             'FixedDt':Variable('flag',0.,{},False,False,'solver_params','to set time step to ref_dt'),
             'ref_dt':Variable('parameter',1.,{},False,False,'solver_params','reference time step'),
             'lattice': Variable('parameter','D2Q9',{},False,False,'solver_params','lattice type for the model'),
             'f': Variable('dist_func',0,{},False,True,'domain_params','LB distribution function'),
             'time': Variable('parameter',0,{},False,True,'solver_params','simulation time in phyiscal units'),
             'iters': Variable('parameter',0,{},False,True,'solver_params','number of iterations in simulations'),
            }
            
    def __init__(self,domain,domain_params,bc_params,solver_params):
        """
        Initializes Incompressible Navier-Stokes equation intance
        """
        domain=self._convert_domaintype(domain)
        if ('rhof' in domain_params) and ('rho' not in domain_params):
            domain_params['rho']=domain_params['rhof']
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
        self._tau = self.visc2tau(visc=self._visc)
            
    def _get_conv_factors(self,*args):
        """
        computes conversion factors for different variables
        """
        if len(args)==3:
            domain,domain_params,solver_params=args
        klist=['rhof','visc','tauref','VelBased','VelRatio','FixedDt','ref_dt','es2']
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
            dx = getattr(self,dx)
        rhof,visc,tauref,VelBased,VelRatio,FixedDt,ref_dt,es2 = vals
        visc_ref = np.max(visc)
        convfactors={}
        if VelBased:
            dt= VelRatio/dx
        elif FixedDt:
            dt=ref_dt
        else:
            visclb=self.tau2visc(tau = tauref)
            visc0= visc_ref/visclb
            dt=dx**2/visc0
        base = {}
        base['L']= dx
        base['M']=rhof*base['L']**3
        base['T']= dt
        convfactors = deepcopy(base)
        for k,v in self._vars.items():
            if v.isphyvar:
                cfact = 1
                for p,v in v.dimension.items():
                    cfact*=base.get(p,1)**v  
                convfactors[k]=cfact
        self.convfactors= convfactors
                
    def tau2visc(self,**kwargs):
        """
        Gets kinematic viscosity (in LB units) from tau
        
        Parameters
        ----------
        tau: int, float or ndarray
            relaxation parameter
        
        Returns
        -------
        visc:    int, float or ndarray
            viscosity in LB units
        """
        if 'tau' in kwargs:
            tau = kwargs['tau']
        else:
            tau = self.tau
        es2=getattr(self,'es2',self._vars['es2'].default)
        return es2*(tau-0.5)
        
    def visc2tau(self,**kwargs):
        """
        gets tau from viscosity (in LB units)

        Parameters
        ----------
        visc: int,float or ndarray
            in LB units

        Returns
        -------
        tau:    int, float or ndarray
            relaxation parameter             
        """
        if 'visc' in kwargs:
            visc = kwargs['visc']
        else:
            visc = self._visc
        es2=getattr(self,'es2',self._vars['es2'].default)
        return visc/es2+0.5
        
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
                    assert lattice == 'D2Q9'
                except AssertionError:
                    raise ValueError("Only D2Q9 lattice implementation exists for 2D")
        elif int(lattice[1]) == 3 and d == 2:
            self._vars['d'].default=d=3
            try:
                assert lattice == 'D3Q19'
                self._vars['q'].default=19
            except AssertionError:
                raise ValueError("Only D3Q19 lattice implementation exists for 3D")
        elif int(lattice[1]) == 3 and d == 3:
            try:
                assert lattice == 'D3Q19'
                self._vars['q'].default=19
            except AssertionError:
                raise ValueError("Only D3Q19 lattice implementation exists for 3D")
        elif int(lattice[1]) == 2 and d == 3:
            self._vars['lattice'].default=lattice='D3Q19'
            self._vars['q'].default=19
        if d==2:
            ns = getattr(_solvers,'ns2d')
        elif d==3:
            ns = getattr(_solvers,'ns3d')   
        #check for forcing models
        if 'forcing_model' in solver_params:
            solver_params['forcing_model'] = solver_params['forcing_model'].lower()
            if collision_model == 'srt':
                av_forcing_schemes =['mlga','sc','guo']
            elif collision_model == 'trt':
                av_forcing_schemes = ['mlga','guo']
            if solver_params['forcing_model'] not in av_forcing_schemes:
                warnings.warn("specified %s forcing model not available set to default 'guo' forcing"%solver_params['forcing_model'])
                solver_params['forcing_model'] = 'guo'
        #compute_macro_var
        self.compute_macro_var = Fmethod(self,ns.compute_macro_var,
                                        ['_f', '_rho0', '_rho', '_prel', '_u','nodetype'])
        #compute edf
        self.compute_edf = Fmethod(self,ns.get_feq,
                                        ['nodetype','_rho', '_u'],inplace_update=False)
        #collison
        if collision_model == 'srt':
            self.collide = Fmethod(self,ns.collide_srt,
                                   ['_f', '_rho',  '_u','_Fv','nodetype','_tau','_forcing_model'])
        elif collision_model == 'trt' :
            self.collide = Fmethod(self,ns.collide_trt,
                                   ['_f', '_rho',  '_u','_Fv','nodetype','_tau','_magic_para','_forcing_model'])        
        else:
            raise NotImplementedError
        #stream
        self.stream = Fmethod(self,ns.stream_and_bounce,['_f', 'nodetype','_periodicity'])      
        #apply_bc
        self.apply_bc = Fmethod(self,ns.apply_bc,
        ['_f', '_rho0','nodetype','_interp'],['_bc'])
    
    @property
    def p(self):
        """
        absolute pressure
        """
        return self.p0 + self.prel
        
   
    def advance(self):
        self.time +=self.dt
        self.iters+=1
        self.collide()
        self.stream()
        self.apply_bc()
        self.compute_macro_var()
    
    @property
    def dt(self):
        """
        """
        return self.convfactors['T']
    
    @property
    def u_mag(self):
        if self.d == 2:
            return (self.u[0]**2 + self.u[1]**2)**0.5
        if self.d == 3:
            return (self.u[0]**2 + self.u[1]**2+ self.u[2]**2)**0.5
    
    def steady_state(self,verbose=True,tol =1e-6,check_iters= 100,err_cal_method='l2norm',max_iters = 50000,skip_boundaries=True):
        while 1:
            rho_u_old = self._rho * self._u
            if skip_boundaries:
                if self.d==2:
                    rho_u_old = rho_u_old[:,1:-1,1:-1]
                elif self.d==3:
                    rho_u_old = rho_u_old[:,1:-1,1:-1,1:-1]  
            self.run(iters = self.iters+check_iters)
            rho_u  = self._rho * self._u
            if skip_boundaries:
                if self.d==2:
                    rho_u = rho_u[:,1:-1,1:-1]
                elif self.d==3:
                    rho_u = rho_u[:,1:-1,1:-1,1:-1]  
            err=getattr(self,err_cal_method)(rho_u, rho_u_old)
            if (err/tol) <1: 
                break
            if self.iters > max_iters:
                warnings.warn("solution did not converged maximum of %s iterations reached"%max_iters)
                break
            if verbose:
                print('Iteration No:%s, Error: %s'%(self.iters,err))
                
        
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
            name=self._bc[i+'_type']  
            if name == 'p':
                val=deepcopy(self._bc[i+'_val'])
                val[0]+= self.p0
                val[0]*=self.convfactors['prel']
            elif name == 'u':
                val=deepcopy(self._bc[i+'_val'])
                val=list(np.array(val)/self.convfactors['u'])
            elif name=='laminar_flow':
                val=deepcopy(self._bc[i+'_val'])
                val[0]/=self.convfactors['u']
            else:
               val=self._bc[i+'_val']               
            bcparams[i]=[name,val]
        return bcparams
        
    @bc.setter
    def bc(self,bc_params):
        """
        setter for bc
        """
        if self.d==2:
            bc=defaultdict(lambda:['nothing',[0,0]])
        else:
            bc=defaultdict(lambda:['nothing',[0,0,0]])
        if hasattr(self,'_bc'):
            bc.update(self.bc)
        bc.update(deepcopy(bc_params))
        if self.d == 2:
            dirs = ['left','right','top','bottom']
        elif self.d == 3:
            dirs = ['left','right','top','bottom','front','back']
        self._bc={}
        for name,val in bc.items():
            if (len(val)==1) and (name in dirs):
                if self.d ==2:
                    val.append([0,0])
                else:
                    val.append([0,0,0])
            elif (len(val)==2):
                if type(val[1])=='float' or type(val[1])=='int':
                    if self.d == 2:
                        val[1]=[val[1],0]
                    elif self.d == 3:
                        val[1]=[val[1],0,0]         
                if (type(val[1]).__name__=='list') or (type(val[1]).__name__=='tuple'):
                    if len(val[1])== 1 and self.d==2:
                        val[1]=[val[1][0],0]
                    elif len(val[1])<3 and self.d==3:
                        if len(val[1])==2:val[1]=[val[1][0],val[1][1],0]
                        if len(val[1])==1:val[1]=[val[1][0],0,0]
        for i in dirs:
            bc[i][1]=list(bc[i][1])
            self._bc[i+'_type']=bc[i][0]  
            if self._bc[i+'_type']=='p':
                bc[i][1][0]-= self.p0
                bc[i][1][0]/=self.convfactors['prel']
                self._bc[i+'_val']=bc[i][1]   
            if self._bc[i+'_type']=='rho':
                bc[i][1][0]/=self.convfactors['rho']
                self._bc[i+'_val']=bc[i][1]       
            elif self._bc[i+'_type']=='u':
                self._bc[i+'_val']=list(np.array(bc[i][1])/self.convfactors['u']) 
            elif self._bc[i+'_type']=='laminar_flow':
                bc[i][1][0]/=self.convfactors['u']
                self._bc[i+'_val']=bc[i][1]
            else:
               self._bc[i+'_val']=bc[i][1]                