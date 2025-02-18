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
import os as _os
def _get_mlist():
    mlist = []
    for f in _os.listdir("."):
        if ('.' in f) and ('make' not in f) and (not f.startswith('_')):
            m = f.split('.')[0]
            if m not in mlist: mlist.append(m)
    return mlist

# for _m in _mlist:
#     # __import__(_m, globals(), locals(), [], 1)
#     print("from . import %s"%_m)
from . import ade2d
from . import ade3d
from . import ccd2d
from . import diff2d
from . import diff3d
from . import fdm
from . import multilevel_ade2d
from . import multilevel_ade3d
from . import multilevel_diff2d
from . import multilevel_diff3d
from . import ns2d
from . import ns3d
from . import update2d
from . import update3d
from . import vans2d
from . import vans3d