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
#conversion of yantra from python2 to python3
#=======================================================================================
#import glob
import fnmatch
import os
import subprocess
for root,dirnames,fnames in os.walk('.'):
    for fname in fnmatch.filter(fnames,'*.py'):
        if fname != '2to3.py':
            subprocess.call(['2to3','-w',os.path.join(root,fname)])

for root,dirnames,fnames in os.walk('.'):
    for fname in fnmatch.filter(fnames,'*.py.bak'):
            subprocess.call(['rm',os.path.join(root,fname)])
