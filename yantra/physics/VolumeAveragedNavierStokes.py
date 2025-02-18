#!/usr/bin/python
# -*- coding: utf-8 -*-
#=======================================================================================
#This File is part of Yantra: A lattice Boltzmann method based tool for multiscale/
#multiphyics simulations
#=======================================================================================
#
#Copyright (C) 2016-2019  <Author> Ravi A. Patel <Email> ravee.a.patel@gmail.com
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
from .IncompressibleNavierStokes import IncompressibleNavierStokes
from yantra import __version__,__license__
__name__ = 'VolumeAveragedNavierStokes'
__author__ = 'Ravi Patel'   
__email__ = 'ravee.a.patel@gmail.com'

class VolumeAveragedNavierStokes(IncompressibleNavierStokes, metaclass=LBMeta):
    _signature = 'yantra.physics.VolumeAveragedNavierStokes.VolumeAveragedNavierStokes'
    _vars=deepcopy(IncompressibleNavierStokes._vars)
    _vars['poros']=Variable('scalar',1.,{},False,False,'domain_params','porosity')
    
    @property 
    def u_superficial(self):
        return self.u * self.poros

    @property 
    def u_superficial_mag(self):
        return self.u_mag * self.poros

    @property 
    def u_darcy(self):
        return self.u_superficial

    @property 
    def u_darcy_mag(self):
        return self.u_superficial_mag()
        
    @property 
    def p_avg(self):
        return self.poros * self.p
    
    def advance(self):
        self.time +=self.dt
        self.iters+=1
        self.compute_macro_var()
        self.collide()
        self.stream()
        self.apply_bc()
    
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
            ns = getattr(_solvers,'vans2d')
        elif d==3:
            ns = getattr(_solvers,'vans3d')   
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
                                       ['_f', '_rho0', '_rho','_poros', '_prel', '_u','nodetype'])
        #compute edf
        self.compute_edf = Fmethod(self,ns.get_feq,
                                        ['nodetype','_rho','_poros','_u'],inplace_update=False)
        #collison
        if collision_model == 'srt':
            self.collide = Fmethod(self,ns.collide_srt,
                                   ['_f', '_rho','_poros',  '_u','_Fv','nodetype',
                                    '_tau','_forcing_model','_periodicity'])
        elif collision_model == 'trt' :
            self.collide = Fmethod(self,ns.collide_trt,
                                   ['_f', '_rho','_poros','_u','_Fv','nodetype','_tau','_magic_para','_forcing_model','_periodicity'])        
        else:
            raise NotImplementedError
        #stream
        self.stream = Fmethod(self,ns.stream_and_bounce,['_f', 'nodetype','_periodicity'])      
        #apply_bc
        self.apply_bc = Fmethod(self,ns.apply_bc,
        ['_f', '_rho0','_poros','nodetype','_interp'],['_bc'])
