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
from yantra._base import LBMeta,Fmethod,_remove_keys
from .MultilevelAdvectionDiffusion import MultilevelAdvectionDiffusion
from yantra import __version__,__license__
__name__ = 'Diffusion'
__author__ = 'Ravi Patel'   
__email__ = 'ravee.a.patel@gmail.com'

class MultilevelDiffusion(MultilevelAdvectionDiffusion, metaclass=LBMeta):
    """
    Solves advection diffusion equation written as
    ..math::
        \partial_{t}\phic = -\vec{nabla} \phi \frac{\delta}{\zeta} \cdot \vec{j}+ ss
        \vec{j}= - D \vec{\nabla}c 
    where,
    c = concentration ..math::[N^{1}L^{-3}]
    u = velocity ..math::[L^{1}T^{-2}]
    j = flux ..math::[N^{1}L^{-2}T^{-1}]
    D = diffusion coefficient ..math::[L^{2}T^{-1}]
    ss = Source/sink term ..math::[N^{1}L^{-3}T^{-1}]
    
    The bounce back condition in streaming step imposes following condition for all nodetype > 0
    \vec{\nabla}c  \cdot \hat{n} = 0 
    """
    _signature = 'yantra.physics.MultilevelDiffusion.MultilevelDiffusion'
    _vars = _remove_keys(MultilevelAdvectionDiffusion._vars,['u'])
    
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
        if d==2:
            diff = getattr(_solvers,'multilevel_diff2d')
        elif d==3:
            diff = getattr(_solvers,'multilevel_diff3d')   
        #compute_macro_var
        self.compute_macro_var = Fmethod(self,diff.compute_macro_var,
                                        ['_f', '_c', '_flux', '_poros','nodetype','_tau'])
        #compute edf
        self.compute_edf = Fmethod(self,diff.compute_edf,
                                        ['_c', '_cphi','_poros','nodetype'],inplace_update=False)

        #collison
        if collision_model == 'trt':
             self.collide = Fmethod(self,diff.collide,
                                   ['_f', '_c', '_cphi','_poros','nodetype',
                                    '_tau','_magic_para','_ss'])      
        else:
            raise ValueError("Only TRT collision model is available")
        #stream
        self.stream = Fmethod(self,diff.stream_and_bounce,['_f', 'nodetype','_periodicity'])      
        #apply_bc
        self.apply_bc = Fmethod(self,diff.boundary_conditions,
        ['_f','nodetype', '_poros', '_tau', '_interp','grid_type'],['_bc'])
