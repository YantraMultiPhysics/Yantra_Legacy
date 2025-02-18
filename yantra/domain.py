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
#This file contains classes related to initialization of Domain 
#
#=======================================================================================


import yantra
from yantra import visualize as v
from yantra._base import Domain
import numpy as np
__all__=['Domain2D','Domain3D']
__name__='domain'
__version__ = yantra.__version__
__author__ = 'Ravi Patel'   
                
class Domain2D(Domain):
    """ 
    Sets 2D domain for simulation
    
    Attributes
    ----------
    corner: tuple  or list
        co-ordinates of bottom left corner of the domain  
    lengths: tuple or list
        length of domain in x and y directions
    grid_type: str
        either `nodal` or `midway` by default `midway`
    nodetype: ndarray
        represents array that can be used to mark a given node with a number and all the \
        nodes that are then associated with that number can assign same parameters. All the \
        nodetype > 0 are considered inactive nodes in physics and nodetype <= 0 considered as active nodes.
    dx: float
        discretization of grid in the domain
    x: ndarray (1D)
        co-ordinates of nodes in x direction
    y: ndarray (1D)
        co-ordinates of nodes in y direction
    ncells: total number of nodes
    d: int, value equal to 2
        dimension of system
        
    Methods
    -------    
    meshgrid()
        returns co-ordinates of nodes in meshgrid format
    draw_rect(center,lengths,idx=1)
        sets values in nodetype array to idx values for all nodes located inside \
        the rectangle
    draw_circle(center,radius,idx=1)
        sets values in nodetype array to idx values for all nodes loacted inside \
        the rectangle
    mark_interface(self, idx1, idx2, idxval,e=[[ 1, 0, -1,  0],[ 0, 1,  0, -1]])
        marks the interface node of idx1 between idx1 and idx2 as idxval
    visualize()
        Visualization for Domain2D instance
    """
    _signature = 'yantra.domain.Domain2D'
    d = 2
    
    def draw_rect(self, center, lengths, idx=1.):
        """
        sets values in nodetype array to idx values for all nodes located \
        inside the rectangle

        Parameters
        ----------
        center: tuple or list
            center of rectangle
        lengths: tuple or list
            length of rectangle in x and y direction
        idx: float or int, optional (default value is 1)
            value to which the nodetype array should be set to within the rectangle
        
        See also
        --------
        :func:`draw_circle`,`draw_multicoated_circle`
        """
        x,y = self.x, self.y
        lengths_bb = (lengths[0] + self.dx ,lengths[1] + self.dx)
        bb = self._bounding_box(center,lengths_bb)
        xmin, xmax = bb[0]
        ymin, ymax = bb[1]
        xx, yy = self.meshgrid(x[xmin:xmax],y[ymin:ymax])
        marker = 1. * ((xx >= (center[0] - lengths[0]/ 2.)) * 
                       (xx <= (center[0] + lengths[0]/ 2.)) *
                       (yy >= (center[1] - lengths[1]/ 2.)) * 
                       (yy <= (center[1] + lengths[1]/ 2.)))
        self.marker = marker
        self.xx = xx
        self.yy = yy
        self.nodetype[ymin:ymax,xmin:xmax] = idx * (marker > 0) + (
        self.nodetype[ymin:ymax,xmin:xmax] * (marker <= 0))
        return

    def draw_circle(self,center,radius,idx=1.):
        """
        sets values in nodetype array to idx values for all nodes located \
        inside the circle

        Parameters
        ----------
        center: tuple or list
            center of circle
        radius: float,int
            radius of circle
        idx: float or int, optional (default value is 1)
            value to which the nodetype array should be set to within the circle
        
        See also
        --------
        :func:`draw_rect`, `draw_multicoated_circle`
        """
        dia = 2*(radius + self.dx)
        x,y = self.x, self.y
        bb = self._bounding_box(center,(dia,dia))
        xmin, xmax = bb[0]
        ymin, ymax = bb[1]
        xx, yy = self.meshgrid(x[xmin:xmax],y[ymin:ymax])
        marker = 1. * ((xx-center[0])**2+ (yy-center[1])**2<=radius**2)
        self.nodetype[ymin:ymax,xmin:xmax] = idx * (marker > 0) + (
        self.nodetype[ymin:ymax,xmin:xmax] * (marker <= 0))
        return

    def draw_multicoated_circle(self,center,radiuses,idxes):
        """
        sets values in nodetype array to idxes values for all nodes located \
        inside the respective layer of multi-coated circle

        Parameters
        ----------
        center: tuple or list
            center of multi-coated circle
        radiuses: tuple, list or ndarray (1D)
            list of radius of multi-coated circles
        idxes: tuple, list or ndarray (1D)
            list of idxes for each radius

        See also
        --------
        :func:`draw_rect`, `draw_circle`
        """
        try:
            assert (len(radiuses)==len(idxes))
        except AssertionError:
            ValueError("both radiuses and idxes should be of same length")
        for radius,idx in zip(radiuses,idxes):
            self.draw_circle(center,radius,idx)
        
    def mark_interface(self, idx1, idx2, idxval,e=[[ 1, 0, -1,  0],[ 0, 1,  0, -1]]):
        """
        marks the interface node of idx1 between idx1 and idx2 as idxval
        
        Parameters
        ----------
        idx1: float
            idx voxel which has to be marked
        idx2: float
            interface node of idx1 between idx1 and idx2 is marked
        idxval: float
            The interface node is marked with idxval
        e: nested list with length equal to 2, optional (default D2Q4 lattice)
            neighbourhood lattice direction zeroth dimension represents x direction \
            first dimension represents y direction
        See also
        --------
        :func:`draw_circle`,`draw_rect`,`draw_multicoated_circle`
        """
        try:
            assert(len(e)==2)
            assert(len(e[0])==len(e[1])>=4)
        except AssertionError:
            ValueError("length of e should be 2 and both list in e should be \
            of same length and contain atleast 4 directions")
        fbb=0
        for ex,ey in zip(e[0],e[1]):
            temp = 1.*(self.nodetype==idx2)
            temp = np.roll(temp,ex,axis=1)
            temp = np.roll(temp,ey,axis=0)
            fbb+= temp
        fbb=(fbb>0)
        cdn = fbb * (self.nodetype==idx1)
        self.nodetype = self.nodetype *(cdn<=0) + idxval*(cdn>0) 


    def draw_ellipse(self,center,radius,aspect_ratio,rotation_angle=0,idx=1.):
        """
        sets values in nodetype array to idx values for all nodes located \
        inside the circle

        Parameters
        ----------
        center: tuple or list
            center of circle
        radius: float,int
            radius of ellipse in  x direction
        aspect_ratio: float,int
            ratio of radius of ellipse in y direction to x direction
        rotation_angle: float,int
            rotation angle for ellipse in radians
        idx: float or int, optional (default value is 1)
            value to which the nodetype array should be set to within the circle
        
        See also
        --------
        :func:`draw_rect`, `draw_multicoated_circle`
        """
        lmax = 2*(max(radius, radius * aspect_ratio) + self.dx)
        x,y = self.x, self.y
        bb = self._bounding_box(center,(lmax,lmax))
        xmin, xmax = bb[0]
        ymin, ymax = bb[1]
        xx, yy = self.meshgrid(x[xmin:xmax],y[ymin:ymax])
        xp = xx-center[0]
        yp = yy-center[1]
        cos = np.cos
        sin = np.sin
        theta =rotation_angle
        xR = xp*cos(theta)   + sin(theta)* yp  
        yR = yp*cos(theta)   - sin(theta)* xp
        x1 = radius
        x2 = radius * aspect_ratio
        marker = 1. * (((xR/x1)**2+ (yR/x2)**2)<=1)
        self.nodetype[ymin:ymax,xmin:xmax] = idx * (marker > 0) + (
        self.nodetype[ymin:ymax,xmin:xmax] * (marker <= 0))
        return
    
    def draw_generalized_ellipse(self,center,radius,aspect_ratio,rotation_angle=0,sw=2,idx=1.):
        """
        sets values in nodetype array to idx values for all nodes located \
        inside the circle

        Parameters
        ----------
        center: tuple or list
            center of circle
        radius: float,int
            radius of ellipse in  x direction
        aspect_ratio: float,int
            ratio of radius of ellipse in y direction to x direction
        rotation_angle: float,int
            rotation angle for ellipse in radians
        sw: float,int
            power of the generalized ellipsoid
        idx: float or int, optional (default value is 1)
            value to which the nodetype array should be set to within the circle
        
        See also
        --------
        :func:`draw_rect`, `draw_multicoated_circle`
        """
        lmax = 2*(max(radius, radius * aspect_ratio) + self.dx)
        x,y = self.x, self.y
        bb = self._bounding_box(center,(lmax,lmax))
        xmin, xmax = bb[0]
        ymin, ymax = bb[1]
        xx, yy = self.meshgrid(x[xmin:xmax],y[ymin:ymax])
        xp = xx-center[0]
        yp = yy-center[1]
        cos = np.cos
        sin = np.sin
        theta =rotation_angle
        xR = xp*cos(theta)   + sin(theta)* yp  
        yR = yp*cos(theta)   - sin(theta)* xp
        x1 = radius
        x2 = radius * aspect_ratio
        marker = 1. * (((np.abs(xR/x1))**sw+ (np.abs(yR/x2))**sw)<=1)
        self.nodetype[ymin:ymax,xmin:xmax] = idx * (marker > 0) + (
        self.nodetype[ymin:ymax,xmin:xmax] * (marker <= 0))
        return

    def visualize(self):
        """
        Visualization for Domain2D instance
        """
        v.visualize_domain(self)
        
class Domain3D(Domain):
    """ 
    Sets 3D domain for simulation
    
    Attributes
    ----------
    corner: tuple  or list
        co-ordinates of bottom left back corner of the domain  
    lengths: tuple or list
        length of domain in x and y directions
    grid_type: str
        either `nodal` or `midway` by default `midway`
    nodetype: ndarray
        represents array that can be used to mark a given node with a number and all the \
        nodes that are then associated with that number can assign same parameters. All the \
        nodetype > 0 are considered inactive nodes in physics and nodetype <= 0 considered as active nodes.
    dx: float
        discretization of grid in the domain
    x: ndarray (1D)
        co-ordinates of nodes in x direction
    y: ndarray (1D)
        co-ordinates of nodes in y direction
    z: ndarray (1D)
        co-ordinates of nodes in z direction
    ncells: total number of nodes
    d: int, value equal to 3
        dimension of system
        
    Methods
    -------    
    meshgrid()
        returns co-ordinates of nodes in meshgrid format
    draw_rect(center,lengths,idx=1)
        sets values in nodetype array to idx values for all nodes located inside \
        the rectangle
    draw_circle(center,radius,idx=1)
        sets values in nodetype array to idx values for all nodes loacted inside \
        the rectangle
    mark_interface(idx1,idx2)  
        sets....
    """
    _signature = 'yantra.domain.Domain3D'
    d = 3
    
    def draw_box(self,center,lengths,idx=1.):
        """
        sets values in nodetype array to idx values for all nodes located \
        inside the box

        Parameters
        ----------
        center: tuple or list
            center of box
        lengths: tuple or list
            length of box in x, y and z direction
        idx: float or int, optional (default value is 1)
            value to which the nodetype array should be set to within the box
        
        See also
        --------
        :func:`draw_sphere`
        """
        x,y,z = self.x, self.y, self.z
        lengths_bb = (lengths[0] + self.dx ,lengths[1] + self.dx,lengths[2] + self.dx)
        bb = self._bounding_box(center,lengths_bb)
        xmin, xmax = bb[0]
        ymin, ymax = bb[1]
        zmin, zmax = bb[2]
        xx, yy, zz = self.meshgrid(x[xmin:xmax],y[ymin:ymax],z[zmin:zmax])
        marker = 1. * ((xx >= (center[0] - lengths[0]/ 2.)) * 
                       (xx <= (center[0] + lengths[0]/ 2.)) *
                       (yy >= (center[1] - lengths[1]/ 2.)) * 
                       (yy <= (center[1] + lengths[1]/ 2.)) *
                       (zz >= (center[2] - lengths[2]/ 2.)) * 
                       (zz <= (center[2] + lengths[2]/ 2.)))
        self.nodetype[zmin:zmax,ymin:ymax,xmin:xmax] =(idx * (marker > 0) +
        self.nodetype[zmin:zmax,ymin:ymax,xmin:xmax] * (marker <= 0))
        return

    def draw_sphere(self, center, radius, idx=1):
        """
        sets values in nodetype array to idx values for all nodes located \
        inside the sphere

        Parameters
        ----------
        center: tuple or list
            center of sphere
        radius: float,int
            radius of sphere
        idx: float or int, optional (default value is 1)
            value to which the nodetype array should be set to within the sphere
        
        See also
        --------
        :func:`draw_box`
        """
        x,y,z = self.x, self.y,self.z
        dia = 2*(radius + self.dx)
        bb = self._bounding_box(center,(dia,dia,dia))
        xmin, xmax = bb[0]
        ymin, ymax = bb[1]
        zmin, zmax = bb[2]
        xx, yy,zz = self.meshgrid(x[xmin:xmax],y[ymin:ymax],z[zmin:zmax])
        marker = 1. * ((xx-center[0])**2+ (yy-center[1])**2+(zz-center[2])**2<=radius**2)
        self.nodetype[zmin:zmax,ymin:ymax,xmin:xmax] = idx * (marker > 0) + (
        self.nodetype[zmin:zmax,ymin:ymax,xmin:xmax] * (marker <= 0))
        return
    
    def draw_ellipsoid(self,center,radius,aspect_ratio=[1,1],rotation_angle=(0,0,0),idx=1):
        """
        """
        x,y,z = self.x, self.y,self.z
        lmax = 2*(max(radius, radius * aspect_ratio[0],radius * aspect_ratio[1]) + self.dx/2)
        bb = self._bounding_box(center,(lmax,lmax,lmax))
        xmin, xmax = bb[0]
        ymin, ymax = bb[1]
        zmin, zmax = bb[2]
        xx, yy,zz = self.meshgrid(x[xmin:xmax],y[ymin:ymax],z[zmin:zmax])
        alpha,beta, gamma = rotation_angle[0],rotation_angle[1],rotation_angle[2]
        #translate
        xp = xx-center[0]
        yp = yy-center[1]
        zp = zz-center[2]
        cos = np.cos
        sin = np.sin
        #rotate (see rotation matrix wikipedia)
        Rx= np.array([[1,         0,          0],
                       [0,cos(alpha),-sin(alpha)],
                       [0,sin(alpha), cos(alpha)]])
        Ry= np.array([[cos(beta) ,0,sin(beta)],
                       [         0,1,        0],
                       [-sin(beta),0,cos(beta)]])
        Rz= np.array([[cos(gamma),-sin(gamma),0],
                       [sin(gamma), cos(gamma),0],
                       [         0,          0, 1]])
        R = Rz.dot(Ry).dot(Rx)
        xr = R[0,0]*xp + R[0,1] *yp + R[0,2] * zp 
        yr = R[1,0]*xp + R[1,1] *yp + R[1,2] * zp 
        zr = R[2,0]*xp + R[2,1] *yp + R[2,2] * zp 
        x1 = radius
        x2 = radius * aspect_ratio[0]
        x3 = radius * aspect_ratio[1]
        marker = 1. * (((xr/x1)**2+(yr/x2)**2+ (zr/x3)**2)<=1)
        self.nodetype[zmin:zmax,ymin:ymax,xmin:xmax] = idx * (marker > 0) + (
        self.nodetype[zmin:zmax,ymin:ymax,xmin:xmax] * (marker <= 0))
        return 

    def draw_generalized_ellipsoid(self,center,radius,aspect_ratio=[1,1],rotation_angle=(0,0,0),sw=2,idx=1):
        """
        """
        x,y,z = self.x, self.y,self.z
        lmax = 2*(max(radius, radius * aspect_ratio[0],radius * aspect_ratio[1]) + self.dx/2)
        bb = self._bounding_box(center,(lmax,lmax,lmax))
        xmin, xmax = bb[0]
        ymin, ymax = bb[1]
        zmin, zmax = bb[2]
        xx, yy,zz = self.meshgrid(x[xmin:xmax],y[ymin:ymax],z[zmin:zmax])
        alpha,beta, gamma = rotation_angle[0],rotation_angle[1],rotation_angle[2]
        #translate
        xp = xx-center[0]
        yp = yy-center[1]
        zp = zz-center[2]
        cos = np.cos
        sin = np.sin
        #rotate (see rotation matrix wikipedia)
        Rx= np.array([[1,         0,          0],
                       [0,cos(alpha),-sin(alpha)],
                       [0,sin(alpha), cos(alpha)]])
        Ry= np.array([[cos(beta) ,0,sin(beta)],
                       [         0,1,        0],
                       [-sin(beta),0,cos(beta)]])
        Rz= np.array([[cos(gamma),-sin(gamma),0],
                       [sin(gamma), cos(gamma),0],
                       [         0,          0, 1]])
        R = Rz.dot(Ry).dot(Rx)
        xr = R[0,0]*xp + R[0,1] *yp + R[0,2] * zp 
        yr = R[1,0]*xp + R[1,1] *yp + R[1,2] * zp 
        zr = R[2,0]*xp + R[2,1] *yp + R[2,2] * zp 
        x1 = radius
        x2 = radius * aspect_ratio[0]
        x3 = radius * aspect_ratio[1]
        marker = 1. * (((np.abs(xr/x1))**sw+(np.abs(yr/x2))**sw+ (np.abs(zr/x3))**sw)<=1)
        self.nodetype[zmin:zmax,ymin:ymax,xmin:xmax] = idx * (marker > 0) + (
        self.nodetype[zmin:zmax,ymin:ymax,xmin:xmax] * (marker <= 0))
        return 
                
    def mark_interface(self, idx1, idx2,idxval,e=[[ 1, 0, -1,  0, 0,  0],
                          [ 0, 1,  0, -1, 0,  0],[ 0, 0,  0,  0, 1, -1]]):
        """
        marks the interface node of idx1 between idx1 and idx2 as idxval
        
        Parameters
        ----------
        idx1: float
            idx voxel which has to be marked
        idx2: float
            interface node of idx1 between idx1 and idx2 is marked
        idxval: float
            The interface node is marked with idxval
        e: nested list with length equal to 3, optional (default D3Q6 lattice)
            neighbourhood lattice direction zeroth dimension represents x direction \
            first dimension represents y direction and second dimension represents z direction
        """
        try:
            assert(len(e)==3)
            assert(len(e[0])==len(e[1])==len(e[2])>=6)
        except AssertionError:
            ValueError("length of e should be 3 and both list in e should be of \
            same length and contain atleast 6 directions")
        fbb=0
        for ex,ey,ez in zip(e[0],e[1],e[2]):
            temp = 1.*(self.nodetype==idx2)
            temp = np.roll(temp,ex,axis=2)
            temp = np.roll(temp,ey,axis=1)
            temp = np.roll(temp,ez,axis=0)
            fbb+= temp
        fbb=(fbb>0)
        cdn = fbb * (self.nodetype==idx1)
        self.nodetype = self.nodetype *(cdn<=0) + idxval*(cdn>0) 

    def visualize(self):
        """
        Visualization for Domain2D instance
        """
        v.visualize_domain(self)