# -*- coding: utf-8 -*-
__doc__="""
The demonstration of restarting AdvectionDiffusion. This demo can also be used to 
test diffusion, multilevelAdvectionDiffusion and multilevelDiffusion
"""
#%% import modules
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
#%%create yantra domain instance
domain = yantra.Domain3D((0,0,0),(200,20,20),1, grid_type = 'nodal')
#%% Define physics
cb = 1
ux=0.01
D=0.1
#domain params 
domain_params={}
domain_params['u']= (ux,0,0)
domain_params['D'] = D
#set boundary conditions
#concentration gradient along x-axis all other boundary periodic
bc_params ={}
bc_params['left']=['c' ,cb]
bc_params['right']=['c',0.0]
#solver parameters
solver_params={}
solver_params['lattice']='D3Q7'
solver_params['collision_model']='trt' #other options 'trt' and 'diff_vel'
solver_params['Dref']= domain_params['D']/10
solver_params['magic_para']= 1./6.
#create physics instance
ade = yantra.AdvectionDiffusion(domain,domain_params,bc_params,solver_params)
#%%run model
ade.run(iters=100)
plt.figure()   
plt.imshow(ade.c[10,:,:])
#%%save model and delete everything
yantra.save(ade,'ade_test')
del ade
del domain
#%%load model
ade=yantra.load('ade_test.ymp')
plt.figure()
plt.imshow(ade.c[10,:,:])
#%%run loaded model 
ade.run(iters=200) # next 100 iterations
plt.figure()
plt.imshow(ade.c[10,:,:])
os.remove('ade_test.ymp')
  
