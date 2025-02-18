
#%% import modules 
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import numpy as np
import random
#%% generation algorithm 
l=50
rmax=5
sw=20
aspect_ratio= [1,0.1]
domain = yantra.Domain3D((0,0,0),(l,l,l),0.25, grid_type = 'nodal',periodicity={'x':1,'y':1,'z':1})
vfrac= 0 
#seed
#for i in range(10):
#    theta=[0,0,0]
#    theta[0]  = random.uniform(0,2*np.pi)
#    theta[1]  = random.uniform(0,2*np.pi)
#    theta[2]  = random.uniform(0,2*np.pi)
#    x = random.uniform(0,l)
#    y = 0
#    z = random.uniform(0,l)
#    domain.draw_generalized_ellipsoid((x,y,z),rmax,aspect_ratio=aspect_ratio,rotation_angle=theta,sw=sw)

#for i in range(10):
#    theta=[0,0,0]
#    theta[0]  = random.uniform(0,2*np.pi)
#    theta[1]  = random.uniform(0,2*np.pi)
#    theta[2]  = random.uniform(0,2*np.pi)
#    x = random.uniform(0,l)
#    y = l
#    z = random.uniform(0,l)
#    domain.draw_generalized_ellipsoid((x,y,z),rmax,aspect_ratio=aspect_ratio,rotation_angle=theta,sw=sw)
#
for i in range(30):
    theta=[0,0,0]
    theta[0]  = random.uniform(0,2*np.pi)
    theta[1]  = random.uniform(0,2*np.pi)
    theta[2]  = random.uniform(0,2*np.pi)
    x = random.uniform(0,l)
    y = random.uniform(0,l)
    z = random.uniform(0,l)
    domain.draw_generalized_ellipsoid((x,y,z),rmax,aspect_ratio=aspect_ratio,rotation_angle=theta,sw=sw)

while vfrac<0.6:
    theta=[0,0,0]
    theta[0]  = random.uniform(0,2*np.pi)
    theta[1]  = random.uniform(0,2*np.pi)
    theta[2]  = random.uniform(0,2*np.pi)
    while 1:
        x = random.uniform(0,l)
        y = random.uniform(0,l)
        z = random.uniform(0,l)
        k = np.argmin(abs(domain.x - x))
        j = np.argmin(abs(domain.y - y))
        i = np.argmin(abs(domain.z - z))
#        if random.uniform(0,1)<=0.05: break
        if domain.nodetype[i,j,k]==1: break
    periodic=  False
    xp = x
    yp = y
    zp = z
    if (x + rmax) > l:
        xp = x-l
        periodic= True
    elif (x- rmax) <0:
        xp = x + l
        periodic= True
    if (y + rmax) > l:
        yp = y-l
        periodic= True
    elif (y - rmax) <0:
        yp = y + l
        periodic= True        
    if (z + rmax) > l:
        zp = z-l
        periodic= True
    elif (z - rmax) <0:
        zp = z + l
        periodic= True 
    domain.draw_generalized_ellipsoid((x,y,z),rmax,aspect_ratio=aspect_ratio,rotation_angle=theta,sw=sw)
    if periodic:
        domain.draw_generalized_ellipsoid((xp,yp,zp),rmax,aspect_ratio=aspect_ratio,rotation_angle=theta,sw=sw)
    vfrac=np.mean(domain.nodetype>0)
    print(vfrac)
yantra._hl.imageToVTK('test_sheets',cellData={'nodetype':domain.nodetype})