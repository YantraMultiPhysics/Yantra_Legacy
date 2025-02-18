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

import numpy as np
from copy import deepcopy
from scipy import ndimage
import yantra
from collections import defaultdict
from yantra import __version__,__license__
#import new
__all__=['Domain', 'LB', 'Variable','Fmethod','LBMeta', 'delvar', 'getvar', 'setvar']
__name__ = '_base'
__author__ = 'Ravi Patel'   
__email__ = 'ravee.a.patel@gmail.com'

def _remove_keys(d,klist):
    """
    Removes keys present in klist from dictonary d
    
    Parameters
    ----------
    d: dict
        dictonary to be operated on
    klist: list
        list of keys to be removed
    
    Returns
    -------
    dict
        d without keys in klist 
    """
    d = deepcopy(d)
    for k in list(d.keys()):
        if k in klist: del d[k]
    return d

def setvar(attr):
    """
    method to set property for physical variables
    """
    def set_var(self,val):
        if self._vars[attr[1:]].isphyvar:
            cfact = self.convfactors[attr[1:]]
        else:
            cfact = 1
        if self._vars[attr[1:]].type=='vector':
            if (type(val).__name__ == 'tuple') or  (type(val).__name__ == 'list'):
                ones = np.ones(self.scalar_shape)
                try:
                    assert(len(val)== self.d)
                except AssertionError:
                    ValueError("%s is a vector and should contain %s as zeroth dimension"
                    %(attr,val))
                val =  list(val)
                for i in range(self.d):
                    val[i]*=ones
                val = np.asfortranarray(np.array(val))
            elif (type(val).__name__ == 'ndarray'):
                try:
                    assert(val.shape== self.vector_shape)
                except AssertionError:
                    ValueError("%s is a vector and should contain %s as zeroth dimension"
                    %(attr,val))
            else:
                ones = np.ones(self.vector_shape)
                val = val * ones               
            val = np.asfortranarray(val/cfact)               
        elif self._vars[attr[1:]].type=='scalar':
            ones = np.ones(self.scalar_shape)
            val = val * ones/cfact   
            val = np.asfortranarray(val) 
        elif self._vars[attr[1:]].type=='dist_func':
            if val.shape == self.dist_func_shape:
                val=np.asfortranarray(val) 
            else:
                ones = np.ones(self.dist_func_shape)
                val = val * ones         
                val = np.asfortranarray(val)   
        elif self._vars[attr[1:]].type=='flag':
            val=bool(val)
        elif self._vars[attr[1:]].type=='int':
            val=int(val)
        else:
            unsupported=['list','str','dict']
            if type(val).__name__ not in unsupported:
                val = val/cfact
        setattr(self,attr,deepcopy(val))
    return set_var

def getvar(attr):
    """
    method to set property for physical variables
    """
    def get_var(self):
        if self._vars[attr[1:]].isphyvar:
            return getattr(self,attr)*self.convfactors[attr[1:]]
        else:
            return deepcopy(getattr(self,attr))          
    return get_var

def delvar(attr):
    """
    method to correctly delete physical variables
    """
    def del_var(self):
        delattr(self,attr)
        return
    return del_var

class Blank(object): pass

class LBMeta(type):
    def __init__(cls, name, bases, dct):
        """
        meta class of the LB created to initalize physical variable as property
        """
        for k,v in cls._vars.items():
            setattr(cls, k,property(fset=setvar("_"+k), fget=getvar("_"+k),
                             fdel=delvar("_"+k)))

class  Variable:
    _args_name = ['type','default','dimension','isphyvar','skip_first_time','belongs_to','doc']
    def __init__(self,*args):
        for name,val in zip(self._args_name,args):
            setattr(self,name,val)
            
        
class Fmethod(object):
    """
    A wrapper to link fortran methods to class methods of objects
    """
    def __init__(self,instance,method,arglist,kwarglist={},outattrs=[],inplace_update=True):
        """
        A wrapper to link fortran methods to class methods of objects
        
        Parameters
        ----------
        instance: any instance 
            instance of any class to which fortran method has to be added
        method: fortran method
            fortran method to be added
        arglist: list
            list of arguments for fortran method
        outattrs: list Default None
            list of output attributes of the instance to be updated by the fortran method     
        """
        self.instance = instance
        self.method = method
        self.arglist = arglist
        self.kwarglist = kwarglist
        self.outattrs = outattrs
        self.inplace_update = inplace_update
        
    def get_args(self):
        """
        gets list of arguments for instance needed for fortran method
        """
        args=[]
        for arg in self.arglist:
            args.append(getattr(self.instance,arg))
        return args

    def get_kwargs(self):
        """
        """
        kwargs={}
        for kwarg in self.kwarglist:
            if type(getattr(self.instance,kwarg)).__name__ == 'dict':
                kwargs.update(getattr(self.instance,kwarg))
            else:
                kwargs.update({kwarg:getattr(self.instance,kwarg)})
        return kwargs
        
    def __call__(self):
        """
        calls the fortran method
        """
        args=self.get_args()
        kwargs=self.get_kwargs()
        if len(self.outattrs) ==0:
            return self.method(*args,**kwargs)
        else:
            out = self.method(*args,**kwargs)
            if self.inplace_update:
                for var,val in zip(self.outattrs,out):
                    setattr(self.instance,var,out)
            else:
                return out

class ExtMethod(Fmethod):
    pass
           
class Domain(object):
    """
    Generic parent class for creating a simulation domain 
    """
    _signature = 'yantra._base.Domain'
    def __init__(self,corner,lengths,dx,grid_type='midway',nodetype=None,periodicity=defaultdict(lambda:0)):
        """
        Initialises Domain instance
        
        Parameters
        ----------
        corner: tuple or list
            corner of the domain for 2D its left bottom corner for 3D it is back left bottom corner
        lengths: tuple or list
            length of domain in x,y directions for 2D and x,y,z directions in 3D
        dx: float
            discretization in physical units
        grid_type: str, optional
             can be specified as `midway` or `nodal`. By default it is set to `midway`
             In `midway` nodes are shifted by half the discretization so that each
             node represents volume around it and boundary nodes are located outside
             the domain of interest thus locating the boundary of the domain  in between boundary node
             and first interior node. This is analogous to cell-centered finite volume meshing.
             In `nodal` the first node represents boundary node and start or end of the domain. Nodes
             are not assumed to represent a volume.             
        nodetype: ndarray, optional
            represents array that can be used to mark a given node with a number and all the
            nodes that are then associated with that number can assign same parameters. All the
            nodetype > 0 are considered inactive nodes in physics and nodetype <= 0 considered as active nodes.
        periodicity: dict, optional
            periodicity is flag if true than domain is treated periodic in that direction
        """
        #check input
        try:
            assert(len(corner)== self.d)
            assert(len(lengths)== self.d)
        except AssertionError:
            ValueError("Corner or lengths argument should contain %s values"%self.d)
        self.corner = corner
        self.lengths = lengths
        self.periodicity=periodicity
        try:
            assert(grid_type.lower()=='nodal' or grid_type.lower()=='midway')
        except AssertionError:
            ValueError("grid_type can either be nodal or midway, %s given"%grid_type)        
        self.grid_type = grid_type.lower()
        self.dx = dx
        self._process_input()
        if np.all(nodetype == None):
            self.nodetype = nodetype
        if nodetype is None:
            self.nodetype = -1. * np.ones(self.shape, order='F')
        else:
            self.nodetype = np.asfortranarray(nodetype)
        
    def _process_input(self):
         """
         processes input for domain
         """
         d,dx,corner,lengths,grid_type = self.d,self.dx,self.corner,self.lengths,self.grid_type
         cord = []
         for x,l in zip(corner,lengths):
             if grid_type == 'nodal':
                 cord.append(np.linspace(x, l, round((l-x)/dx)+1))
             elif grid_type == 'midway':
                 x= x+dx/2; l = l-dx/2
                 cord.append(np.linspace(x, l, round((l-x)/dx)+1))
         if d == 2:
             self.x, self.y = cord
             self.y = self.y[::-1]
             self.shape = (len(self.y),len(self.x))
             self.ny, self.nx = self.shape
         elif d == 3:
             self.x,self.y,self.z = cord
             self.y = self.y[::-1]
             self.z = self.z[::-1]
             self.shape = (len(self.z),len(self.y),len(self.x))
             self.nz,self.ny,self.nx = self.shape
             
    def meshgrid(self,*args):
         """
         creates co-ordinates in meshgrid format
         
         Parameters
         ----------
         x,y: list or 1D ndarray (for 2D), optional
         x,y,z: list or 1D ndarray (for 3D), optional
         
         Returns
         -------
         x,y: ndarray (for 2D)
         x,y,z: ndarray (for 3D)
         """
         if self.d == 2:
             if len(args) > 0 and len(args)<=2:
                 try:
                     x,y  = args
                 except IndexError:
                     ValueError("There should be two input lists")
             else:
                 x,y = self.x,self.y
             xx,yy = np.meshgrid(x,y)
             return np.asfortranarray(xx),np.asfortranarray(yy)
         elif self.d == 3:
             if len(args) > 0 and len(args)<=3:
                 try:
                     x,y,z  = args
                 except IndexError:
                     ValueError("There should be three input lists")
             else:
                 x,y,z = self.x,self.y,self.z
             yy,zz,xx = np.meshgrid(y,z,x)
             return np.asfortranarray(xx),np.asfortranarray(yy),np.asfortranarray(zz)
    

    def percolated_clusters(self,idx,stencil_name,axis=None):
        """
        returns ndarray with value 1 for percolated cluster of phase equal to idx \
        in nodetype and 0 for non-percolated cluster of phase equal to idx in \ 
        nodetype
        
        parameters
        ----------
        idx: int,float
            value of idx in nodetype for which percolated clusters is needed
        stencil_name:str
            name of stencil to be considered for getting percolated cluster
            for 2D: 'D2Q5' and 'D2Q9'
            for 3D: 'D3Q7','D3Q19' and 'D3Q27'
        axis: str, optional (None)
            axis direction along with percolated clusters are evaluated. If none \
            1 is assigned to cluster percolated in any boundary direction
        
        returns
        -------
        ndarray
            ndarray equal to 1 for percolated clusters of phase equal to idx in nodetype \
            and zero for non-percolated cluster
        
        see also
        --------
        :func:`degree_percolation`
        """
        d=self.d
        nodetype = self.nodetype
        stencil = self._get_stencil_map(stencil_name)
        lblarray, nlabel = ndimage.label(nodetype == idx, structure=stencil)
        #x axis
        if d==2:
            nlabel = lblarray[:, 0] * (lblarray[:, 0] == lblarray[:, -1])
        elif d==3:
            nlabel = lblarray[:, :, 0] * (lblarray[:, :, 0] == lblarray[:, :, -1])
        nlabel = np.unique(nlabel)
        pclusters_x = np.zeros(np.shape(lblarray))
        for i in nlabel:
            if i > 0: pclusters_x += 1. * (lblarray == i)
        #y axis
        lblarray, nlabel = ndimage.label(nodetype == idx, structure=stencil)  
        if d==2:
            nlabel = lblarray[0, :] * (lblarray[0, :] == lblarray[ -1, :]) 

        elif d==3:
            nlabel = lblarray[:, 0, :] * (lblarray[:, 0, :] == lblarray[:, -1, :])
        nlabel = np.unique(nlabel)
        pclusters_y = np.zeros(np.shape(lblarray))
        for i in nlabel:
            if i > 0: pclusters_y += 1. * (lblarray == i)
        #zaxis
        if d == 3:
            nlabel = lblarray[0, :, :] * (lblarray[0, :, :] == lblarray[-1, :, :])
            nlabel = np.unique(nlabel)
            pclusters_z = np.zeros(np.shape(lblarray))
            for i in nlabel:
                if i > 0: pclusters_z += 1. * (lblarray == i)
        del nlabel, lblarray
        if axis=='x': return 1.*(pclusters_x>0)
        if axis=='y': return 1.*(pclusters_y>0) 
        if axis=='z': return 1.*(pclusters_z>0)
        if axis==None and d==2: return 1.*((pclusters_x+pclusters_y)>0)
        if axis==None and d==3: return 1.*((pclusters_x+pclusters_y+pclusters_z)>0)

    def labelled_clusters(self,idx,stencil_name,axis=None):
        """
        returns ndarray with labelled clusters and array containing label values
        
        parameters
        ----------
        idx: int,float
            value of idx in nodetype for which percolated clusters is needed
        stencil_name:str
            name of stencil to be considered for getting percolated cluster
            for 2D: 'D2Q5' and 'D2Q9'
            for 3D: 'D3Q7','D3Q19' and 'D3Q27'
        axis: str, optional (None)
            axis direction along with percolated clusters are evaluated. If none \
            1 is assigned to cluster percolated in any boundary direction
        
        returns
        -------
        ndarray
            ndarray with values equal to label number  
        ndarray
            list of labels of clusters
        see also
        --------
        :func:`degree_percolation`
        """
        nodetype = self.nodetype
        stencil = self._get_stencil_map(stencil_name)
        lblarray, nlabel = ndimage.label(nodetype == idx, structure=stencil)
        return lblarray, nlabel
    
    def degree_percolation(self,idx,stencil_name,axis=None):
        """
        degree of percolation of phase equal to idx in nodetype
        
        parameters
        ----------
         idx: int,float
            value of idx in nodetype for which percolated clusters is needed
        stencil_name:str
            name of stencil to be considered for getting percolated cluster
            for 2D: 'D2Q5' and 'D2Q9'
            for 3D: 'D3Q7','D3Q19' and 'D3Q27'
        axis: str, optional (None)
            axis direction along with percolated clusters are evaluated. If none \
            1 is assigned to cluster percolated in any boundary direction
        
        returns
        -------
        float
            degree of percolation of phase equal to idx in nodetype
        """
        vsum = np.sum(self.nodetype==idx)
        percolated = self.percolated_clusters(idx,stencil_name,axis)
        psum = np.sum(percolated)
        return psum/vsum
    
    def _bounding_box(self,center,lengths):
        """get bounding box"""
        d,dx = self.d,self.dx
        if d == 2:
            cord = [self.x,self.y]
        elif d == 3:
            cord = [self.x,self.y,self.z]
        bb = []
        for i in range(d):
            mini=np.abs(cord[i]-((center[i]-lengths[i]/2.)-dx)).argmin()
            maxi=np.abs(cord[i]-((center[i]+lengths[i]/2.)+dx)).argmin() 
            if i==0:
                bb.append((mini,maxi+1))
            else:
                bb.append((maxi,mini+1))
        return bb
        
    @staticmethod
    def _get_stencil_map(name):
        """get stencil for percolation"""
        if name == 'D3Q27':
            return [[[1, 1, 1],
                     [1, 1, 1],
                     [1, 1, 1]],
                    [[1, 1, 1],
                     [1, 1, 1],
                     [1, 1, 1]],
                    [[1, 1, 1],
                     [1, 1, 1],
                     [1, 1, 1]]]
        if name == 'D3Q7':
            return [[[0, 0, 0],
                     [0, 1, 0],
                     [0, 0, 0]],
                    [[0, 1, 0],
                     [1, 1, 1],
                     [0, 1, 0]],
                    [[0, 0, 0],
                     [0, 1, 0],
                     [0, 0, 0]]]
        if name == 'D3Q19':
            return [[[0, 1, 0],
                     [1, 1, 1],
                     [0, 1, 0]],
                    [[1, 1, 1],
                     [1, 1, 1],
                     [1, 1, 1]],
                    [[0, 1, 0],
                     [1, 1, 1],
                     [0, 1, 0]]]
        if name == 'D2Q5':
            return [[0, 1, 0],
                    [1, 1, 1],
                    [0, 1, 0]]
        if name == 'D2Q9':
            return [[1, 1, 1],
                    [1, 1, 1],
                    [1, 1, 1]]
    @property
    def ncells(self):
        """
        Total number of cells in domain
        """
        if self.d  == 2:
            return self.nx * self.ny
        if self.d == 3:
            return self.nx * self.ny * self.nz
    @classmethod
    def use_img(cls,img,dx=1,periodicity=defaultdict(lambda:0)):
       corner = (0,0,0)
       lengths = list(img.shape)
       lengths.reverse()
       lengths = [i*dx for i in lengths]
       return cls(corner,lengths,dx,'midway',img,periodicity)
    
    @classmethod
    def loaddict(cls,toload):
        """
        loads the domain data from a dictonary generated using savedict
        """
        return cls(toload['corner'],toload['lengths'],toload['dx'],toload['grid_type'],
                   toload['nodetype'],toload['periodicity'])
    
    def savedict(self):
        """
        saves data needed to load domain using loaddict
        """
        return {'corner':self.corner,
                'lengths':self.lengths,
                'dx':self.dx,
                'grid_type':self.grid_type,
                'nodetype':self.nodetype,
                'periodicity':self.periodicity,
                'signature': self._signature
                }

    def whoami(self):
        """
        returns signature of class
        """
        print(self._signature)
        
class LB(object, metaclass=LBMeta):
    """
    Generic parent class for lattice Boltzmann model implementation
    """
    _vars={}
    _signature = 'yantra._base.LB'
    
    def  __init__(self, domain,domain_params={},bc_params={},solver_params={}):
        """
        Initalizes LB object
        """
        domain=self._convert_domaintype(domain)
        self._get_conv_factors(domain,domain_params,solver_params)
        self._read_input(domain,domain_params,bc_params,solver_params)

    def _get_conv_factors(self, domain,domain_params,solver_params):
        """
        Computes conversion factors and stores in convfactors
        """
        raise NotImplementedError    

    @classmethod
    def loaddict(cls,toload):
        """
        loads the model from dictionary generated using savedict
        """
        return cls(toload['domain'],toload['domain_params'],
                   toload['bc_params'],toload['solver_params'])
    
    def savedict(self):
        """
        saves information of model generated necessary to load it using loaddict
        """
        #get domain data
        vars_in_domain = ['nx','ny','nz','d','grid_type','nodetype','dx','x','y','z','periodicity']        
        domain = {}
        for var in vars_in_domain:
            try:
                domain[var]=getattr(self,var)
            except AttributeError:
                pass
        #get domain_params and solver parameter data
        domain_params={}
        solver_params={}
        for k,v in self._vars.items():
            if v.belongs_to == 'domain_params':
                domain_params[k]=getattr(self,k)
            elif v.belongs_to == 'solver_params':
                solver_params[k]=getattr(self,k)
        tosave = {'domain':domain,'domain_params':domain_params,
                 'solver_params':solver_params,'bc_params':self.bc,
                 'signature':self._signature}
        return tosave        
    
    @staticmethod
    def _convert_domaintype(domain):
        """
        quick hack to convert domain saved in dict for in savedict method to a class as 
        its difficult to pickle a class
        """
        vars_from_domain = ['nx','ny','nz','d','grid_type','nodetype','dx','x','y','z','periodicity']        
        if type(domain).__name__ == 'dict':
            temp = Blank()
            for var in vars_from_domain:
                try:
                    setattr(temp,var,domain[var])
                except:KeyError
            domain = temp
        return domain
        
    def _read_input(self, domain,domain_params,bc_params,solver_params):
        """
        reads input and sets relevant attributes to the instance 
        """
        #read domain
        vars_from_domain = ['nx','ny','nz','grid_type','nodetype','dx','x','y','z','periodicity']        
        for var in vars_from_domain:
            try:
               val= getattr(domain,var)
               setattr(self,var,val)
            except: NameError
        #copy meshgrid method of Domain class
#        meshgrid = object.__getattribute__(Domain,'meshgrid')
#        self.meshgrid = new.instancemethod(meshgrid,self,self.__class__)
        self.d = solver_params.get('d',self._vars['d'].default)
        self.q = solver_params.get('q',self._vars['q'].default)
        try:
            assert (int(domain.d) == int(self.d))
        except AssertionError:
            ValueError(
                        """
                        Dimensions of domain and physics should be same.  
                        Got d=%s for domain and d=%s for physics. 
                        """%(domain.d,self.d)
                        )
        self.periodicity = domain.periodicity
        if self.d == 2: 
            self.scalar_shape  = (self.ny, self.nx)
            self.vector_shape = (self.d,self.ny,self.nx)
            self.dist_func_shape = (self.ny, self.nx, self.q)
        elif self.d == 3:
            self.scalar_shape = (self.nz, self.ny, self.nx)
            self.vector_shape = (self.d, self.nz, self.ny, self.nx)
            self.dist_func_shape = (self.nz, self.ny, self.nx, self.q)
        if domain.grid_type == 'nodal':
            self._vars['interp'].default = 0
        elif domain.grid_type == 'midway':
            self._vars['interp'].default = 1           
        #read domain_params and solver_params
        inputs = deepcopy(domain_params); inputs.update(solver_params)
        for k,v in self._vars.items():
            if v.skip_first_time == False or  v.skip_first_time == None:
                if callable(v.default):
                    setattr(self,k,inputs.get(k,v.default()))
                else:
                    setattr(self,k,inputs.get(k,v.default))
        #assign bc
        self.bc = bc_params
        #initialize remaining variables skiped during first time 
        for k,v in self._vars.items():
             if v.skip_first_time == True:
                if callable(v.default):
                    setattr(self,k,inputs.get(k,v.default()))
                else:
                    setattr(self,k,inputs.get(k,v.default))

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
    def ncells(self):
        """
        Total number of cells in domain
        """
        if self.d  == 2:
            return self.nx * self.ny
        if self.d == 3:
            return self.nx * self.ny * self.nz
            
    @property
    def periodicity(self):
        if self.d == 2:
            return {'x':self._periodicity[0],'y':self._periodicity[1]}
        if self.d == 3:
            return {'x':self._periodicity[0],'y':self._periodicity[1],'z':self._periodicity[1]}
    
    @periodicity.setter
    def periodicity(self,p):
        if self.d == 2:
            self._periodicity = [p['x'],p['y']]
            return
        elif self.d==3:
            self._periodicity = [p['x'],p['y'],p['z']]
            return
        
    def meshgrid(self,*args):
         """
         creates co-ordinates in meshgrid format
         
         Parameters
         ----------
         x,y: list or 1D ndarray (for 2D), optional
         x,y,z: list or 1D ndarray (for 3D), optional
         
         Returns
         -------
         x,y: ndarray (for 2D)
         x,y,z: ndarray (for 3D)
         """
         if self.d == 2:
             if len(args) > 0 and len(args)<=2:
                 try:
                     x,y  = args
                 except IndexError:
                     ValueError("There should be two input lists")
             else:
                 x,y = self.x,self.y
             xx,yy = np.meshgrid(x,y)
             return np.asfortranarray(xx),np.asfortranarray(yy)
         elif self.d == 3:
             if len(args) > 0 and len(args)<=3:
                 try:
                     x,y,z  = args
                 except IndexError:
                     ValueError("There should be three input lists")
             else:
                 x,y,z = self.x,self.y,self.z
             yy,zz,xx = np.meshgrid(y,z,x)
             return np.asfortranarray(xx),np.asfortranarray(yy),np.asfortranarray(zz)
         
    def whoami(self):
        """
        returns signature of class
        """
        print(self._signature)

    @staticmethod
    def l1norm(a,b):
        return (np.sum(np.abs(a-b))/np.sum(np.abs(a)))
    
    @staticmethod
    def l2norm(a,b):
        return (np.sum((a-b)**2)/np.sum(a**2))**0.5

    @staticmethod
    def max_abs_diff(a,b):
        return np.nanmax(np.abs(a-b))
    
class Multicomponent(object):
    """
    base class for multicomponent transport models
    """
    _signature = 'yantra._base.Multicomponent'
    def __init__(self,eqn,components,domain,domain_params,bc_params,solver_params):
        """
        """
        #initialize instances of eqn
        self.components = components 
        self.eqn = eqn 
        eqn = getattr(yantra,eqn)
        for name in components:
            my_domain_params = deepcopy(domain_params)
            my_bc_params = deepcopy(bc_params)
            my_solver_params = deepcopy(solver_params)  
            if name in my_domain_params:
                my_domain_params.update(my_domain_params[name])
            if name in my_bc_params:
                my_bc_params.update(my_bc_params[name])
            if name in my_solver_params:
                my_solver_params.update(my_solver_params[name])
            for _name in components:
                if (_name in my_domain_params) and (_name!=name):
                    del my_domain_params[_name]
                if (_name in my_solver_params) and (_name!=name):
                    del my_solver_params[_name]
                if (_name in my_bc_params) and (_name!=name):
                    del my_bc_params[_name]
            setattr(self,name,eqn(domain,my_domain_params,my_bc_params,my_solver_params))
    
    def call(self,method,*args,**kwargs):
        """
        calls a method for all components
        """
        output={}
        for name in self.components:
            model = getattr(self,name)
            f= getattr(model,method)
            output[name]=f(*args,**kwargs)
        return output

    def set_attr(self,param,val,component_dict=True):
        for name in self.components:
            model = getattr(self,name)
            if component_dict:
                setattr(model,param,val[name])
            else:
                setattr(model,param,deepcopy(val))

    def get_attr(self,param):
        out={}
        for name in self.components:
            model = getattr(self,name)
            out[name]= deepcopy(getattr(model,param))
        return out

    @classmethod
    def loaddict(cls,toload):
        pass
    
    def savedict(self):
        pass

class NewMultiComponent(dict):
    """
    base class for multicomponent transport models
    """
    _signature = 'yantra._base.Multicomponent'
    def __init__(self,*args):
        if len(args)>1:
            self._init_new_multicomponent(*args)
        else:
            self._load_multicomponent(*args)
        
    def _load_multicomponent(self,toload):
        self.components = toload['components']
        self.eqn = toload['eqn']
        eqn = getattr(yantra,self.eqn)
        for name in self.components:
            self[name]= eqn.loaddict(toload[name])
    
    def _init_new_multicomponent(self,eqn,components,domain,domain_params,bc_params,solver_params):
        """
        """
        #initialize instances of eqn
        self.components = components
        self.eqn = eqn 
        eqn = getattr(yantra,eqn)
        for name in components:
            my_domain_params = deepcopy(domain_params)
            my_bc_params = deepcopy(bc_params)
            my_solver_params = deepcopy(solver_params)  
            if name in my_domain_params:
                my_domain_params.update(my_domain_params[name])
            if name in my_bc_params:
                my_bc_params.update(my_bc_params[name])
            if name in my_solver_params:
                my_solver_params.update(my_solver_params[name])
            for _name in components:
                if (_name in my_domain_params) and (_name!=name):
                    del my_domain_params[_name]
                if (_name in my_solver_params) and (_name!=name):
                    del my_solver_params[_name]
                if (_name in my_bc_params) and (_name!=name):
                    del my_bc_params[_name]
            self[name]=eqn(domain,my_domain_params,my_bc_params,my_solver_params)
    
    def call_comp_method_all(self,method,*args,**kwargs):
        """
        calls a method for all components
        """
        output={}
        for name in self.components:
            f= getattr(self[name],method)
            output[name]=f(*args,**kwargs)
        return output
        
    def call_comp_method(self,name,method,*args,**kwargs):
        """
        call a method for a given component
        """
        f = getattr(self[name],method)
        return f(*args,**kwargs)
        
    def set_comp_attr(self,name,param,val):
        setattr(self[name],param,deepcopy(val))
    
    def set_comp_attr_all(self,param,val,same_values=False):
        for name in self.components:
            if same_values:
                setattr(self[name],param,deepcopy(val))
            else:
                setattr(self[name],param,deepcopy(val[name]))

    def get_comp_attr(self,name,param):
        return deepcopy(getattr(self[name],param))

    def get_comp_attr_all(self,param):
        out={}
        for name in self.components:
            out[name]= deepcopy(getattr(self[name],param))
        return out
        
    @classmethod
    def loaddict(cls,toload):
        return cls(toload)
    
    def savedict(self):
        tosave = {}        
        tosave['components']=self.components
        tosave['eqn']=self.eqn
        for name in self.components:
            tosave[name]=self[name].savedict()
        return tosave 