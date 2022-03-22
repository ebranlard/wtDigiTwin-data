""" 
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
create=False
# create=True
runSim=True

tMax = 10
tMax = None

# modelName = 'B111111_hydro0'; fstFilename = 'Spar_MD0HD1SD0_F111111/Main.fst'; hydro=True

# modelName = 'B111111_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F111111/Main.fst'; hydro=True


modelName = 'B101010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F101010/Main.fst'; hydro=True
# modelName = 'B001000_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F001000/Main.fst'; hydro=True
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010/Main.fst'; hydro=True 
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_F000010_NoHydro/Main.fst'; hydro=True

# modelName = 'B101010_hydro0'; fstFilename = 'MD0HD1SD0_F101010/Main.fst'; hydro=True

# Hydro False (i.e. see StructSim..)
#modelName = 'B100010_hydro0'; fstFilename = 'MD0HD0SD0_F100010/Main.fst'; hydro=False
#modelName = 'B100010'; fstFilename = 'MD0HD0SD0_F100010/Main.fst'; hydro=False
#modelName = 'F100010T0RNA_fnd_noLoads'; fstFilename = 'MD0HD0SD0_F100010/Main.fst'; hydro=False


# --- Generate python package
if create:
    if modelName[0]=='B':
        generateOneRigidBodyModel(modelName)
    else:
        generateModel(modelName, aero_forces=False, moor_loads=False, hydro_loads=True)

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

    sim.setCoupledHydroLoads()
    #uop=None
    #du=None
    qop  = None
    qdop = None
    # Using mean as op
    # qop  = np.array([np.mean(dfFS[c]) for c in WT.q_channels])
    # qdop = np.array([np.mean(dfFS[c]) for c in WT.qd_channels])*0
    # --- Simulation
    #sim.setInputs(u, du, uop, qop, qdop)
    sim.simulate(out=True, prefix='_hydroPyCoupled')
    sim.plot(export=True, prefix='_hydroPyCoupled')

    # --- Plot forcing
    #fig,axes = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    #axes=np.atleast_1d(axes)
    #sim.sysNL.plot_forcing(fig=fig, axes=axes, includeCK=False, label='Non linear', c='k')
    #sim.sysLI.plot_forcing(fig=fig, axes=axes, includeCK=True, plotCK0=True, label='Linear')
    #axes[0].legend()
plt.show()
