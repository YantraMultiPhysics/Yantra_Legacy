
#%% import modules 
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import random
#%% generation algorithm 
l=50
rmax=2
sw=2
aspect_ratio= [1,10]
x=l/2
y=l/2
z=l/2
domain = yantra.Domain3D((0,0,0),(l,l,l),0.25, grid_type = 'nodal',periodicity={'x':1,'y':1,'z':1})
theta=[np.pi/4.,0,0]
domain.draw_generalized_ellipsoid((x,y,z),rmax,aspect_ratio=aspect_ratio,rotation_angle=theta,sw=sw)
yantra._hl.imageToVTK('template',cellData={'nodetype':domain.nodetype})