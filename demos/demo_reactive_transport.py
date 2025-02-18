"""
Reactive transport modeling using Yantra-xGems
"""
#%%

import sys,os
myxgemspath = '/home/ravi/ownCloud/codes/xgems/build/lib' #change this to path where xgems is compiled. If directory structure is same as in repo dont change this
sys.path.append(myxgemspath)
from xgems import ChemicalEngine
import numpy as np
from scipy.special import erfc
import matplotlib.pylab as plt
from copy import deepcopy as copy
import time
#set matplotlib font size and type
font = {'family' : 'serif',
         'size'   : 34}
plt.rc('font',**font)
#%%
class VarDef(object):
    """
    a class to define variable used for convenience 
    """
    _args_name=['type','default','doc']
    def __init__(self,*args):
        for name,val in zip(self._args_name,args):
            setattr(self,name,val)
#%%
class xGems(object):
    """
    xGems object for coupling with reactive transport
    """
    _vars={'gems_init_file':VarDef('str','','Gems file for initialization ends with .lst'),
           'gems_ic_files':VarDef('list',[],'Gems dbr files list for setting initial conditions'),
           'ic_cell_labels':VarDef('scalar',0,'an array specify cell label corresponding to initial state given by dbr file'),
           'gems_bc_files':VarDef('list',[],'Gems dbr files for setting boundary conditions'),
           'T':VarDef('scalar',0,'temprature for each cell'),
           'pe':VarDef('scalar',0,'pe'),
           'P':VarDef('scalar',0,'pressure for each cell'),
           'ncells':VarDef('param',0,'number of cells'),
            '_poros':VarDef('scalar',1,'internal variable for porosity'),
            '_phaseVolFrac':VarDef('dict',1,'internal for  phase volume fractions'),
            '_phaseConc':VarDef('dict',1,'internal for  phase concentration'),
            'pH':VarDef('scalar',0,'pH'),
            'pE':VarDef('scalar',0,'pe'),
            }
    
    def __init__(self,ncells,shape,inputs={}):
        """
        initalizes gems reaction module
        
        Input
        -----
         inputs: dict
             dictionary containing inputs for initializing xGems classs
             inputs['gems_init_file']: input file obtained from GEMS3K to initialize xGems
             inputs['gems_ic_files']: list containing files obtained from GEMS3K for initial conditions
             inputs['gems_bc_files']:  list containing files obtained from GEMS3K for boundary conditions
             inputs['ic_cell_labels']: array to specify label for cell corresponding to index of files in initial conditions
        """
        self.ncells = ncells  
        self.shape = shape
        self._read_inputs(inputs)  
        self.gem = ChemicalEngine(self.gems_init_file) 
        self.nelements = self.gem.numElements()
        self.element_names = []
        for i in range(self.nelements):
            self.element_names.append(self.gem.elementName(i))
        self.nphases= self.gem.numPhases()
        self.phase_names =[]
        for i in range(self.nphases):
            self.phase_names.append(self.gem.phaseName(i))
        
    def _read_inputs(self,inputs):
        """
        a convenience method to read inputs
        """
        skip_items=['ncells']
        for k,v in self._vars.items():
            if k not in skip_items:
                if v.type=='scalar':
                    val = inputs.get(k,v.default)
                    if type(val).__name__ == 'ndarray': val.flatten()
                    setattr(self,k,*np.array([1]*self.ncells))
                else:
                    setattr(self,k,inputs.get(k,v.default))
            
    def advance(self,c_dict):
        """
        advances a timestep in calculation
        
        Input
        -----
        c_dict: dict
            dictionary of concentrations of elements obtained from  transport step
        """
        for name in list(c_dict.keys()):
            self.c_mob[name]=c_dict[name].flatten()
        for i in range(self.ncells):
            T,P = self.T[i],self.P[i]
            b=[]
            for name in self.element_names:
                b.append(c_dict[name][i]+self.c_immob[name][i])
            self.gem.equilibrate(T,P,b)
            T,P,b = self.gem.temperature(),self.gem.pressure(),self.gem.elementAmounts()
            c_mob  = self.gem.elementAmountsInPhase(self.phase_names.index('aq_gen'))
            for j,name in enumerate(self.element_names):
                self.c_mob[name][i]=c_mob[j]
                self.c_immob[name][i]=b[j]-c_mob[j]
            pVol = self.gem.phaseVolumes()
            pAmt  = self.gem.phaseAmounts()
            V = self.gem.systemVolume()
            for j,name in enumerate(self.phase_names):
                self._phVfrac[name][i] = pVol[j]/V
                self._phConc[name][i] = pAmt[j]             
            self._poros[i]=self._phVfrac['aq_gen'][i]
            self.pH[i] = self.gem.pH()
            self.pE[i]=self.gem.pe()
        return 
 
    def get_ic(self):
        """
        provides element concentrations for initializing transport solver 
        
        Returns
        -------
        c_dict: dict
            dictionary of concentrations of elements  
        """
        out=[]
        for fname in self.gems_ic_files:
            T,P,b=self.read_lst_file(fname)
            self.gem.equilibrate(T, P, b)  
            T = self.gem.temperature() 
            P = self.gem.pressure()  
            b = self.gem.elementAmounts()  
            V = self.gem.systemVolume()
            c_mob  = self.gem.elementAmountsInPhase(self.phase_names.index('aq_gen'))
            ph,pe =self.gem.pH(),self.gem.pe() 
            pVol =self.gem.phaseVolumes()
            pAmount = self.gem.phaseAmounts()
            out.append([T,P,V,c_mob,b-c_mob,ph,pe,pVol,pAmount])
        self.c_immob ={}
        self.c_mob = {}
        for name in self.element_names:
            self.c_mob[name]=[]  
            self.c_immob[name]=[]
        self._phConc ={}
        self._phVfrac={}
        self._poros = []
        for name in self.phase_names:
            self._phConc[name]=[]
            self._phVfrac[name]=[]            
        for i in range(self.ncells):
            idx = self.ic_cell_labels[i]
            T,P,V,c_mob,c_immob,ph,pe,phaseVol,phaseAmt = out[idx]
            self.T[i]=T
            self.P[i]=P
            for j,name in enumerate(self.element_names):
                self.c_immob[name].append(c_immob[j])
                self.c_mob[name].append(c_mob[j])
            for j,name in enumerate(self.phase_names):
                self._phConc[name].append(phaseAmt[j])             
                self._phVfrac[name].append(phaseVol[j]/V)
                if name == 'aq_gen':
                    self._poros.append(phaseVol[j]/V)
            self.pH[i] = ph
            self.pE[i]=pe
        return self.get_cdict()
    
    def get_bc(self):
        """
        provides element concentration for boundary conditions 
        """
        out=[]
        idx = 0
        for fname in self.gems_bc_files:
            T,P,b=self.read_lst_file(fname)
            self.gem.equilibrate(T, P, b)
            T = self.gem.temperature() 
            P = self.gem.pressure() 
            b = self.gem.elementAmounts() 
            c_mob  = self.gem.elementAmountsInPhase(self.phase_names.index('aq_gen'))
            out.append({})
            for i,name in enumerate(self.element_names):
                out[idx][name]=c_mob[i]
            idx+=1
        return out

    def get_cdict(self):
        """
        returns element concentration for 
        """
        temp = {}
        for name in self.element_names:
            temp[name]=self.c_mob[name].reshape(self.shape)
        return 

    @property
    def porosity(self):
        """
        porosity
        """
        return self._poros 
   
    @property
    def phase_volume_frac(self):
        """
        phase volume fractions
        """
        return self._phVfrac
    
    @property
    def phase_conc(self):
        """
        phase concentrations
        """
        return self._phConc
    
    @staticmethod
    def read_lst_file(fname):
        gem = ChemicalEngine(fname)
        return gem.temperature(),gem.pressure(),gem.elementAmounts()
#%%        
def test_xgems_reactive_trnasport():
    """
    tests gems reactive transport for benchmark problem
    ref: Shao H et al. Appl. Geochem. 2009 Jul 1, 24(7):1287-300.
    """
    #reactive transport example
    domain=Domain1D(0.5,0.5/80) #set up simulation domain
    inputs={}
    inputs['phi']=0.32
    inputs['u']= 9.375e-6
    inputs['D0'] = 9.375e-6 * 0.0067
    inputs['gems_init_file']=os.path.join('resources','CalciteIC-dat.lst')
    inputs['gems_ic_files']=[os.path.join('resources','CalciteIC-dat.lst')]
    inputs['gems_bc_files']=[os.path.join('resources','CalciteBC-dat.lst')]
    inputs['ic_cell_labels'] = 0
    rxn= xGems(domain,inputs) #initalize xgems
    #get initial and boundary conditions from gems
    ic = rxn.get_ic()
    bc = rxn.get_bc()
    #setup correct_boundary and initial conditions
    for name in rxn.element_names:
        inputs[name]={}
        inputs[name]['c'] = ic[name]
        inputs[name]['bc']= {'left':['c',bc[0][name]],'right':['j',0]}
    trans= MultiComponentTransport(domain,inputs,rxn.element_names) #initialize transport
#    run model we use sequential non-interative approach for coupling 
    i=0
    while trans.t < 21000:
        i+=1
        trans.advance() #run transport 
        cdict = trans.get_cdict() #get conc
        rxn.advance(cdict) #supply conc and run gems
        cdict = rxn.get_cdict() #get new equilibriated concentrations
        trans.update_c(cdict) #update concentration in transpor module
        if i%100==0:print ("Time: %s"%trans.t)
#   post processing of results
    fig,ax1= plt.subplots()
    cdict = trans.get_cdict()
    x = trans.x
    outlist=['Ca','Mg','Cl']
    lineformat = ['-k','-b','-r','-g','-c']
    i=0
    for name in outlist:
        ax1.plot(x,cdict[name],lineformat[i],label = name,linewidth=6)
        i+=1
    ax1.set_ylabel('ionic concentration (M)')
    ax1.set_xlabel('Distance (m)')
    ax1.set_ylim([0,2.1e-3])
    ax1.set_xlim([0,0.5])
    ax2 = ax1.twincells()
    outlist=['Calcite','Dolomite-dis']
    for name in outlist:
        ax2.plot(x,rxn.phase_conc[name],lineformat[i],label = name,linewidth=6)
        i+=1
    ax2.set_ylabel('mineral concentration (M)')
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc=1,fontsize=24)   
    ax2.set_ylim([0,4e-4])
    ax2.set_xlim([0,0.5])
    ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0),useMathText=True)
    ax2.ticklabel_format(style='sci',axis='y',scilimits=(0,0),useMathText=True)
    plt.show()        
    return trans,rxn
    
if __name__ == '__main__':
    trans,rxn=test_xgems_reactive_trnasport()
