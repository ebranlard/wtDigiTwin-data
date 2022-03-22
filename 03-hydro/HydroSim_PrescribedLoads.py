""" 
 

 F000010
	 NOTE: to get the best results, it's best to set the uop of F_hx and M_hy to zero
	       because these variables tend to overshoot at t=0, and the mean is not a representative
	       operating point value
	OR: to have a sufficiently long time series/perturbation about a mean

  F001000
		For improved results set the "OP" of Fz to M g  (since mean might have some offset)
		
		NOTE:  
		  The heave equation might need to be replaced by a simple restoring heave/spring motion
		  Otherwise, things might end up unstable
"""
import os
import numpy as np    
import matplotlib.pyplot as plt
import importlib
# yams
from welib.tools.clean_exceptions import *
from welib.yams.windturbine import FASTWindTurbine
from welib.yams.models.simulator import *
from welib.yams.models.generator_oneRigidBody import generateOneRigidBodyModel
from welib.yams.models.generator import generateModel

# ---- Script parameters
# create=False
create=True
runSim=True

tMax = 10
tMax = None

modelName = 'B101010_hydro0'; fstFilename = 'SparNoRNA_F101010/Main.fst'; hydro=True
modelName = 'B001000_hydro0'; fstFilename = 'SparNoRNA_F001000/Main.fst'; hydro=True
modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_F000010/Main.fst'; hydro=True 
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_F000010_NoHydro/Main.fst'; hydro=True

# modelName = 'B101010_hydro0'; fstFilename = 'MD0HD1SD0_F101010/Main.fst'; hydro=True

# Hydro False (i.e. see StructSim..)
#modelName = 'B100010_hydro0'; fstFilename = 'MD0HD0SD0_F100010/Main.fst'; hydro=False
#modelName = 'B100010'; fstFilename = 'MD0HD0SD0_F100010/Main.fst'; hydro=False
#modelName = 'F100010T0RNA_fnd_noLoads'; fstFilename = 'MD0HD0SD0_F100010/Main.fst'; hydro=False


# --- Generate python package
if create:
    #generateOneRigidBodyModel(modelName)
    generateModel(modelName)

# --- Run non linear and linear simulation using a FAST model as input
if runSim:

    # --- Setup Sim
    print('----------------------- SETUP SIMULATION -----------------------------------------')
    WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)  # TODO
    sim = SimulatorFromOF(WT, modelName=modelName, packageDir='py')
    if modelName[0]=='B':
        time, dfFS, p = sim.setupSim(tMax=tMax, flavor='onebody', J_at_Origin=True)
    else:
        time, dfFS, p = sim.setupSim(tMax=tMax, J_at_Origin=True)

    if not hydro:
        pass # default to 0
    else:
        # --- Inputs
        u=sim.u
        u['F_hx'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFxi_[N]'])
        u['F_hz'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFzi_[N]'])
        u['M_hy'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroMyi_[N-m]'])

        # --- Linear model input operating point
        uop=sim.uop
        uop['F_hx'] = np.mean(dfFS['HydroFxi_[N]'].values)  *0
        uop['F_hz'] = np.mean(dfFS['HydroFzi_[N]'].values)
        uop['F_hz'] = p['M_B']*p['g']
        uop['M_hy'] = np.mean(dfFS['HydroMyi_[N-m]'].values)*0

        # --- Linear pertubation inputs
        du = sim.du
        for i,su in enumerate(sim.info['su']):
            if su=='F_hx': du[i,:] = dfFS['HydroFxi_[N]'].values     - uop['F_hx']
            if su=='F_hz': du[i,:] = dfFS['HydroFzi_[N]'].values     - uop['F_hz']  #- p['M_B']*p['g']
            if su=='M_hy': du[i,:] = dfFS['HydroMyi_[N-m]'].values   - uop['M_hy']

        print('F_hz mean', np.mean(dfFS['HydroFzi_[N]'].values))
        print('M g      ', p['M_B']*p['g'])
    #uop=None
    #du=None
    qop  = None
    qdop = None
    # Using mean as op
    # qop  = np.array([np.mean(dfFS[c]) for c in WT.q_channels])
    # qdop = np.array([np.mean(dfFS[c]) for c in WT.qd_channels])*0
    # --- Simulation
    #sim.setInputs(u, du, uop, qop, qdop)
    sim.simulate()
    sim.plot()

    # --- Plot forcing
    #fig,axes = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    #axes=np.atleast_1d(axes)
    #sim.sysNL.plot_forcing(fig=fig, axes=axes, includeCK=False, label='Non linear', c='k')
    #sim.sysLI.plot_forcing(fig=fig, axes=axes, includeCK=True, plotCK0=True, label='Linear')
    #axes[0].legend()
plt.show()
