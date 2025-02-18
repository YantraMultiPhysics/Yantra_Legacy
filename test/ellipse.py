
#%% import modules 
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import random
#%%create yantra domain instance
domain = yantra.Domain2D((0,0),(100,100),0.05, grid_type = 'midway',periodicity={'x':1,'y':1})
domain.draw_ellipse((0,0),10,0.25,np.pi/4)
#while np.mean(domain.nodetype>0)<=0.3:
#    theta  = random.uniform(0,2*np.pi)
#    x = random.uniform(0,100)
#    y = random.uniform(0,100)
#    domain.ellipse((x,y),8,1,theta)
domain.visualize()
yantra.export_to_vtk(domain,'domain',['nodetype'])
#%%generalized ellipsoid
#domain = yantra.Domain2D((0,0),(500,100),0.5, grid_type = 'nodal',periodicity={'x':1,'y':1})
#while np.mean(domain.nodetype>0)<=0.65:
#    theta  = random.uniform(0,2*np.pi)
#    x = random.uniform(0,100)
#    y = random.uniform(0,100)
#    periodic=  False
#    xp = x
#    yp = y
#    if (x + 8) > 100:
#        xp = x-100
#        periodic= True
#    elif (x - 8) <0:
#        xp = x + 100
#        periodic= True
#    if (y + 8) > 100:
#        yp = y-100
#        periodic= True
#    elif (y - 8) <0:
#        yp = y + 100
#        periodic= True        
#    domain.draw_generalized_ellipse((x,y),2.5,0.5,theta,sw=10)
#    if periodic:
#        domain.draw_generalized_ellipse((xp,yp),2.5,0.5,theta,sw=10)
# 
#domain.visualize()