import matplotlib.pyplot as plt
# from tabulate import tabulate
#%%
solvers = ['ade2d','ade3d','ns2d','ns3d']
nproc=[1,2,4,8,16,24]
plt.style.use('seaborn-darkgrid')
fig=plt.figure(figsize=(10,4))
ax=[]
ax.append(fig.add_subplot(121))
ax.append(fig.add_subplot(122))
ax[0].set_xlabel('number of threads')
ax[0].set_ylabel('speed up')
ax[1].set_xlabel('number of threads')
ax[1].set_ylabel('MLUPS')
#%%
for solver in solvers:
    mlups = []
    for i in nproc:
        fname = solver + '_' + 'np' + '_' + str(i) + '.txt'
        with open(fname,'r') as f:
            val=(f.readline().split(':')[-1])
            mlups.append(float(val))
    mlups_1 = mlups[0]
    # print('-'*4,solver,'-'*4)
    # table = [[nproc[i],mlups[i]] for i in range(len(nproc))]
    # table.insert(0,['processors','mlups'])
    # print(tabulate(table))
    speedup = [i/mlups_1 for i in mlups]
    ax[0].plot(nproc,speedup, '-o',label = solver)
    ax[1].plot(nproc,mlups, '-o',label = solver)
#%%
ax[0].plot(nproc, nproc,label='linear')
ax[0].legend()
#ax[1].legend()
fig.tight_layout()
plt.show()    
fig.savefig('scaling.svg',quality=100, dpi=1000,optimize=True)
