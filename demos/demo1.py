# coding: utf-8
__doc__="""
Creating domain,saving it and loading it...
"""
#%% import modules 
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
##%% 2D domain test
##create a periodic domain
#d = yantra.Domain2D((0,0),(50,50),1,grid_type = 'midway', periodicity = {'x':1,'y':1,'z':1})
##d.draw_circle((25,25),10,1)
##d.draw_circle((0,0),10,1)
##d.draw_circle((0,50),10,1)
##d.draw_circle((50,50),10,1)
##d.draw_circle((50,0),10,1)
#d.draw_rect((25,25),(10,20),1)
#
#plt.figure()
#plt.imshow(d.nodetype)
#yantra.save(d,'domain2d')
#del d
##%% load saved domain
#d=yantra.load('domain2d.ymp')
#plt.figure()
#plt.imshow(d.nodetype)
##%% 3D domain test
##create a periodic domain
d = yantra.Domain3D((0,0,0),(50,50,50),1,grid_type = 'midway', periodicity = {'x':1,'y':1,'z':1})
d.draw_sphere((0,0,0),20,1)
yantra.export_to_vtk(d,'domain',['nodetype'])
#plt.figure()
#plt.imshow(d.nodetype[:,:,25])
#yantra.save(d,'domain3d')
#del d
##%% load saved domain
#d=yantra.load('domain3d.ymp')
#plt.figure()
#plt.imshow(d.nodetype[:,:,25])
#os.remove('domain2d.ymp')
#os.remove('domain3d.ymp')