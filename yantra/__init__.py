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
__author__ = 'Ravi A. Patel'   
__email__ = 'ravee.a.patel@gmail.com'
__license__ = 'GPL V3 and later'
from .version import version as __version__
#import domain,physics
from .domain import *
from .physics import *
import numpy as _np
from ._pyevtk.evtk import hl as _hl
from ._pyevtk.evtk.vtk import VtkGroup
#import _base
import importlib
import pickle as pickle
from scipy import io as  _io
import numpy as _np
import time as _time

def load(fname):
    """
    loads a model from a file
    
    Parameters
    ----------
    fname: str
        
    """
    d=pickle.load(open( fname, 'rb')) 
    return loaddict(d)

def save(m,fname):
    d= savedict(m)
    pickle.dump(d, open( fname+'.ymp', 'wb'), protocol=pickle.HIGHEST_PROTOCOL) 

def loaddict(d):
    """
    loads an yantra's model instance from dictionary
    
    parameters
    ----------
    d: dict
        dict generated using savedict method
    """
    signlist = d['signature'].split('.')
    sign = '.'.join(signlist[:-1])
    m=importlib.import_module(sign)
    m = getattr(m,signlist[-1])
    return m.loaddict(d)

def savedict(m):
    """
    loads an yantra's model instance from dictionary
    
    parameters
    ----------
    d: dict
        dict generated using savedict method
    
    Returns
    -------
    dict
        dictonary contain all information about m to be reloaded
    """
    return m.savedict()    

def export_to_matlab(m,fname,varname):
    """
    saves variables to have access in matlab
    
    Parameters
    ----------
    
    """
    mdict={}
    if m.d == 2:
        mdict['x'],mdict['y']= m.meshgrid()
    elif m.d == 3:
        mdict['x'],mdict['y'],mdict['z']= m.meshgrid()
    for name in varname:
        mdict[name]=getattr(m,name)
    _io.savemat(fname,mdict,long_field_names=True)

def export_to_vtk(module,fname,varname):
    """
    module: Yantra Domain or Physics instance
    fname: file name for input
    varname: list of variables for input
    """
    m  = module
    outvars={}
    if m.d == 2:
        for name in varname:
            val = getattr(m,name)
            if len(val.shape) > m.d: 
                val = list(val)
                val =  [_np.expand_dims(v,0).T for v in val]
                val = [_np.ascontiguousarray(v[:,::-1]) for v in val]
                val.append(val[-1]*0)
                val = tuple(val)
            else:
                val = _np.expand_dims(val,0).T
                val = _np.ascontiguousarray(val[:,::-1])
            outvars[name] = val
        if m.grid_type =='nodal':
            
            _hl.imageToVTK(fname,origin=(m.x[0],m.y[-1],0),spacing=(m.dx,m.dx,0),pointData=outvars)
        else:
            _hl.imageToVTK(fname,origin=(m.x[0],m.y[-1],0),spacing=(m.dx,m.dx,0),cellData=outvars)
    elif m.d == 3:
        for name in varname:
            val = getattr(m,name)
            if len(val.shape) > m.d: 
                val = list(val)
                val = [_np.ascontiguousarray(v.T[:,::-1,::-1]) for v in val]
                val = tuple(val)
            else:
                val = _np.ascontiguousarray(val.T[:,::-1,::-1])                
            outvars[name]=val
        if m.grid_type == 'nodal':
            _hl.imageToVTK(fname,origin=(m.x[0],m.y[-1],m.z[-1]),spacing=(m.dx,m.dx,m.dx),pointData=outvars)
        else:
            _hl.imageToVTK(fname,origin=(m.x[0],m.y[-1],m.z[-1]),spacing=(m.dx,m.dx,m.dx),cellData=outvars)

def dict_to_vtk(fname,origin,spacing,vardict,grid_type='midway'):
    """
    low level functionality to export to vtk 
    """
    if grid_type == 'nodal':
        _hl.imageToVTK(fname,origin=origin,spacing=spacing,pointData=vardict)
    else:
        _hl.imageToVTK(fname,origin=origin,spacing=spacing,cellData=vardict)

        
def l1norm(a,b):
    return (_np.sum(_np.abs(a-b))/_np.sum(_np.abs(a)))

def l2norm(a,b):
    return (_np.sum((a-b)**2)/_np.sum(a**2))**0.5


class timer(object):
    def __init__(self, func = lambda :None):
        """
        Yantra timer class
        """
        self.func = func
        self.t0 = _time.time()
    
    @property
    def current_time(self):
        return _time.time()
    
    @property
    def time_lapsed(self):
        return self.current_time -self.t0
    
    def time_taken(self):
        d,h,m,s = self.to_dhms(self.time_lapsed)
        if (d == 0) and (h==0) and (m==0):
            if s < 1:
                return "Time lapsed (ms): %f"%(s*1e3)
            else:
                return "Time lapsed (s): %f"%(s)
        elif (d == 0) and (h==0):
            return "Time lapsed (mm:ss): %02i:%02i"%(m,s)
        elif (d == 0):
            return "Time lapsed (hh:mm:ss): %02i:%02i:%02i"%(h,m,s)
        else:
            return "Time lapsed (dd:hh:mm:ss): %02i:%02i:%02i:%02i"%(d,h,m,s)

    def to_dhms(self,t):
         m,s = divmod(t,60)
         h,m = divmod(m,60)
         d,h=divmod(h,24)
         return d,h,m,s
     
    def __str__(self):
        return self.time_taken()
    
    def __repr__(self):
        return "Instance for yantra timer class"
    
    def __call__(self,*args,**kwargs):
        res = self.func(*args,**kwargs)
        print(self.time_taken())
        return res