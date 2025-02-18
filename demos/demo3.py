# -*- coding: utf-8 -*-

__doc__="""
The demonstration of restarting IncompressibleNavierStokes.
"""

#%% import modules
import sys,os
PARENT = '..'
sys.path.append(os.path.join(os.path.dirname(__file__), PARENT))
import yantra
import matplotlib.pylab as plt
#%%   setup domain
l = 20
dx=1.0
domain= yantra.Domain2D((0,0),(l,0.5*l),dx, grid_type = 'midway')
#%% set navier stokes
pin =0.001
pout = 0
tauref = 1.0
visc = 1./6.
domain_params={}
domain_params['visc']=visc
domain_params['rho0']=1
bc_params = {}
bc_params['left']=['p',[pin,0]]
bc_params['right']=['p',[pout,0]]
bc_params['top']=['u',[0,0]]
bc_params['bottom']=['u',[0,0]]
solver_params={'lattice':'D2Q9'}
solver_params['tauref']=tauref
ns = yantra.IncompressibleNavierStokes(domain,domain_params,bc_params,solver_params)
#%% run model
ns.run(iters=100)
x,y = ns.meshgrid()
plt.figure()
plt.contourf(x,y,(ns.u[0]**2+ns.u[1]**2)**0.5)
plt.axis('image')
plt.title('velocity vector')
plt.colorbar()
plt.show()
#%%save model
yantra.save(ns,'test_ns')
del domain,ns
#%%load model
#ns = pickle.load(open('test_ns.ymp','rb'))
ns=yantra.load('test_ns.ymp')
plt.figure()
plt.contourf(x,y,(ns.u[0]**2+ns.u[1]**2)**0.5)
plt.axis('image')
plt.title('velocity vector')
plt.colorbar()
plt.show()
#%%run model further
ns.run(iters=1000)
x,y = ns.meshgrid()
plt.figure()
plt.contourf(x,y,(ns.u[0]**2+ns.u[1]**2)**0.5)
plt.axis('image')
plt.title('velocity vector')
plt.colorbar()
plt.show()
os.remove('test_ns.ymp')