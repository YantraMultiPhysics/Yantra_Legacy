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
#This file contains functions to visualize data from different classes 
#
#=======================================================================================

import yantra
__name__='visualize'
__version__ = yantra.__version__
__author__ = 'Ravi Patel'
import numpy as np
import sys
try:
    import matplotlib.pylab as plt
except ImportError:
    _matplotlib = False
try:
    from mayavi import mlab
except ImportError:
    _mlab = False
if sys.version_info > (3,0):
    range = range
else:
    range = xrange
        
def _plot2d(plottype,m,var,title=None,val=None,colorbar=True,mask_solid=False,mask_fluid=False,
            *args,**kwargs):
    """
    """
    x,y = m.meshgrid()
    if np.all(val== None): val = getattr(m,var)
    if mask_solid:  val=np.ma.masked_where(m.nodetype>0,val)
    if mask_fluid:  val=np.ma.masked_where(m.nodetype<=0,val)
    myplot=getattr(plt,plottype)
    myplot(x,y,val,*args,**kwargs)
    if title== None: 
        if getattr(val,'__doc__',None)==None:
            title=var
        else:
            title=val.__doc__
    plt.title(title,fontsize=22)
    plt.tick_params(axis='both', which='major', labelsize=18)
    if colorbar: plt.colorbar().ax.tick_params(labelsize=14)
    plt.axis('image')
    
def _plot3d(plottype,m,var,title=None,val=None,colorbar=True,mask_solid=False,mask_fluid=False,
            *args,**kwargs):
    x,y,z = m.meshgrid()
    extent =[np.min(x),np.max(x),np.min(y),np.max(y),np.min(z),np.max(z)]
    if np.all(val== None): val = getattr(m,var)
    if mask_solid:  
        ntype = list(m.nodetype.flatten())
        x=list(x.flatten());y=list(y.flatten());z=list(z.flatten())
        val = list(val.flatten())
        ext = len(x)
        x = [x[i] for i in range(ext) if ntype[i] <= 0]
        y = [y[i] for i in range(ext) if ntype[i] <= 0]
        z = [z[i] for i in range(ext) if ntype[i] <= 0]
        val = [val[i] for i in range(ext) if ntype[i] <= 0]
    if mask_fluid:  
        ntype = list(m.nodetype.flatten())
        x=list(x.flatten());y=list(y.flatten());z=list(z.flatten())
        val = list(val.flatten())
        ext = len(x)
        x = [x[i] for i in range(ext) if ntype[i] > 0]
        y = [y[i] for i in range(ext) if ntype[i] > 0]
        z = [z[i] for i in range(ext) if ntype[i] > 0]
        val = [val[i] for i in range(ext) if ntype[i] > 0]
    myplot=getattr(mlab,plottype)
    mlab.figure()
    h = mlab.gcf()
    h.scene.background = (1, 1, 1)
    h.scene.foreground = (0, 0, 0)
    myplot(x,y,z,val,*args,**kwargs)
    if colorbar:     mlab.colorbar()
    mlab.outline(figure=h,extent=extent,line_width=3.0,color=(0,0,0))
    mlab.axes(extent =extent)
    mlab.orientation_axes()
    if title== None: 
        if getattr(val,'__doc__',None)==None:
            title=var
        else:
            title=val.__doc__
    mlab.title(title)
    mlab.show()
    return h

def pcolor(m,var,title= None,val=None,colorbar=True,mask_solid=False,mask_fluid=False,*args,**kwargs):
    """
    pcolor plot of the variable
    
    parameters
    ----------
    m: instance of yantra classes
        m can be Domain instance or any instance of classes in Physics module
    var: str
        variable to be plotted
    title: str
        title of the plot
    """
    if m.d == 2:
        figure=plt.figure()
        _plot2d('pcolormesh',m,var,title,val,colorbar,mask_solid,mask_fluid=mask_fluid,*args,**kwargs)
        plt.show() 
        return figure
    elif m.d == 3:
        h= _plot3d('points3d',m,var,title,val,colorbar,mask_solid,mask_fluid=mask_fluid,
                scale_factor=m.dx,mode='cube',scale_mode='scalar',
                line_width=0,*args,**kwargs)
        return h
    
def phaseplot(m,var,phaseid,title=None,cmap='gray',color= (0,0,0)):
    """
    pcolor plot of the variable showing only values equal to phaseid
     
    parameters
    ----------
    m: instance of yantra classes
        m can be Domain instance or any instance of classes in Physics module
    var: str
        variable to be plotted
    phaseid: int,float
        var==phaseid would be shown
    title: str
        title of the plot
    val: optional 
    """
    val = getattr(m,var)
    if title==None:
        if getattr(val,'__doc__',None)==None:
            title=var + ' == %s'%phaseid
        else:
            title=val.__doc__ + ' == %s'%phaseid
    val= 1.*(val==phaseid)
    if m.d == 2:
        figure=pcolor(m,var,title,val,cmap=plt.get_cmap(cmap),
                      colorbar=False,mask_fluid=True)
    elif m.d == 3:
        figure=pcolor(m,var,title,val,color=color,colorbar=False,mask_fluid=True)
    return figure    

def contour(m,var,title=None,isolevels= 40,val=None,mask_solid=False,colorbar=True):
    """
    contour plot of the variable
    
    parameters
    ----------
    m: instance of yantra classes
        m can be Domain instance or any instance of classes in Physics module
    var: str
        variable to be plotted
    title: str
        title of the plot
    isolevels: int, optional (default = 40)
        number of isolevels for contour
    """
    if m.d == 2:
        figure=plt.figure()
        _plot2d('contourf',m,var,title,val,colorbar,mask_solid,isolevels)
        plt.show()
        return figure
    
def  visualize_domain(domain):
    """
    visualization of domain
    
    parameters
    ----------
    domain: Domain2D or Domain3D instance
        instance of domain to be visualized
    """
    if domain.d == 2:        
        figure=pcolor(domain,'nodetype')
        suptitle ='corner:%s lengths: %s, dx: %s'%(domain.corner,domain.lengths,domain.dx)
        figure.suptitle(suptitle,fontsize=22, fontweight='bold')
        plt.show()
    if domain.d == 3:        
        title = "nodetype\ncorner:%s lengths: %s, dx: %s"%(domain.corner,domain.lengths,domain.dx)
        figure=pcolor(domain,'nodetype',title = title)
        title = "nodetype>0\ncorner:%s lengths: %s, dx: %s"%(domain.corner,domain.lengths,domain.dx)
        figure=pcolor(domain,'nodetype',mask_fluid=True,title = title)        
        title = "nodetype<=0\ncorner:%s lengths: %s, dx: %s"%(domain.corner,domain.lengths,domain.dx)
        figure=pcolor(domain,'nodetype',mask_solid=True,title = title)        