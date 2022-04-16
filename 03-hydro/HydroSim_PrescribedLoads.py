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
create=False
# create=True
runSim=True

tMax = 10
tMax = None
# tMax = 64

qop  = None
qdop = None



#---  Hydro False (i.e. see StructSim..)
#modelName = 'B100010_hydro0'; fstFilename = 'MD0HD0SD0_F100010/Main.fst'; hydro=False
#modelName = 'B100010'; fstFilename = 'MD0HD0SD0_F100010/Main.fst'; hydro=False
# modelName = 'F100010T0RNA_fnd_noLoads'; fstFilename = 'MD0HD0SD0_F100010/Main.fst'; hydro=False



# --- Spar No RNA
# - Loads at 0
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010/Main.fst'; hydro=True; # Perfect
# modelName = 'B101010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F101010/Main.fst'; hydro=True
# modelName = 'B001000_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F001000/Main.fst'; hydro=True
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010/Main.fst'; hydro=True 
# modelName = 'B000010_hydroO'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010/Main.fst'; hydro=True 
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_F000010_NoHydro/Main.fst'; hydro=True
# modelName = 'B111111_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F111111/Main.fst'; hydro=True #NL/Lin: perfect except yaw for NL

# - Loads at Origin
# modelName = 'B000010_hydroO'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010/Main.fst'; hydro=True  # NOT OK
# modelName = 'B101010_hydroO'; fstFilename = 'SparNoRNA_MD0HD1SD0_F101010/Main.fst'; hydro=True 
# modelName = 'B111111_hydroO'; fstFilename = 'SparNoRNA_MD0HD1SD0_F111111/Main.fst'; hydro=True #NL/Lin: perfect except yaw for NL (TODO BUT REFz HAS WRONG SIGN)


# --- Spar
# - Loads at 0
# modelName = 'B001000_hydro0'; fstFilename = 'Spar_MD0HD1SD0_F001000/Main.fst'; hydro=True # NL/Lin:Perfect
# modelName = 'B111111_hydro0'; fstFilename = 'Spar_MD0HD1SD0_F111111/Main.fst'; hydro=True # NL/Lin:OK, acc a bit funky, yaw lin off
# - Loads at Origin
# modelName = 'B111111_hydroO'; fstFilename = 'Spar_MD0HD1SD0_F111111/Main.fst'; hydro=True # 


# --- Tetra Spar
# --- Loads at 0
# modelName = 'B001000_hydro0'; fstFilename = 'TS_MD0HD1SD0_F001000/Main.fst'; hydro=True; # ~OK 
# modelName = 'B001000_hydro0'; fstFilename =  'TS_MD0HD1SD0_F001000_NoCdCpCa/Main.fst'; hydro=True # ~OK
# modelName = 'B000010_hydro0'; fstFilename =  'TS_MD0HD1SD0_F000010_NoCdCpCa/Main.fst'; hydro=True ; qop=[-2.266*np.pi/180] 
# modelName = 'B000010_hydroO'; fstFilename =  'TS_MD0HD1SD0_F000010_NoCdCpCa/Main.fst'; hydro=True ; qop=[-2.2*np.pi/180] # NOTE: At O. OK
# modelName = 'B100010_hydro0'; fstFilename = 'TS_MD0HD1SD0_F100010/Main.fst'; hydro=True;  qop=[-1,-2.266*np.pi/180]  # TODO
# modelName = 'B100010_hydro0'; fstFilename = 'TS_MD0HD1SD0_F100010_NoRef/Main.fst'; hydro=True;  qop=[0,-0.48*np.pi/180]  # TODO 
# modelName = 'B101010_hydro0'; fstFilename = 'TS_MD0HD1SD0_F101010/Main.fst'; hydro=True ; qop=[-2.5, 3.2, -2*np.pi/180]  # ~OK, need more mass (2kg)
# modelName = 'B111111_hydro0'; fstFilename = 'TS_MD0HD1SD0_F111111_NoCdCpCa/Main.fst'; hydro=True; # qop=[-2.0,0.08,2.9,0,-0.036,0.0018]; qdop=[0,0,0,0,0,0];
# modelName = 'B111111_hydro0'; fstFilename = 'TS_MD0HD1SD0_F111111_NoCdCpCa_Init/Main.fst'; hydro=True; # qop=[-2.0,0.08,2.9,0,-0.036,0.0018]; qdop=[0,0,0,0,0,0];
modelName = 'B111111_hydro0'; fstFilename = 'TS_MD0HD1SD0_F111111/Main.fst'; hydro=True;  qop=[-2.0,0.08,2.9,0,-2*np.pi/180,0.0*np.pi/180]; qdop=[0,0,0,0,0,0]; # ~OK, except yaw
# --- Loads at Origin
# modelName = 'B101010_hydroO'; fstFilename = 'TS_MD0HD1SD0_F101010/Main.fst'; hydro=True 
# modelName = 'B111111_hydroO'; fstFilename = 'TS_MD0HD1SD0_F111111/Main.fst'; hydro=True; qop=[-2.0,0.08,2.9,0,-2*np.pi/180,0.0*np.pi/180]; qdop=[0,0,0,0,0,0]; # ~OK, except yaw


# --- Generate python package
if create:
    if modelName[0]=='B':
        generateOneRigidBodyModel(modelName, CG_on_z=False)
    else:
        generateModel(modelName, aero_forces=False, moor_loads=False, hydro_loads=True)

# --- Run non linear and linear simulation using a FAST model as input
if runSim:
    # --- Setup Sim
    print('----------------------- SETUP SIMULATION -----------------------------------------')
    WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], algo='OpenFAST')
    sim = SimulatorFromOF(WT, modelName=modelName, packageDir='py')
    if modelName[0]=='B':
        time, dfFS, p = sim.setupSim(tMax=tMax, flavor='onebody', J_at_Origin=True)
    else:
        time, dfFS, p = sim.setupSim(tMax=tMax, J_at_Origin=True)

    if not hydro:
        pass # default to 0
    else:
        sim.setPrescribedHydroInputs(meanMhy=True)
    #uop=None
    #du=None
    sim.qop  = qop
    sim.qdop = qdop
    # Using mean as op
    # sim.qop  = np.array([np.mean(dfFS[c]) for c in WT.q_channels])
    # sim.qdop = np.array([np.mean(dfFS[c]) for c in WT.qd_channels])*0
    # --- Simulation
    #sim.setInputs(u, du, uop, qop, qdop)
    sim.simulate(out=True, prefix='_hydroPyPrescrLoads')
    sim.plot(export=True,  prefix='_hydroPyPrescrLoads', nPlotCols=3)

    # --- Plot forcing
    #fig,axes = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    #axes=np.atleast_1d(axes)
    #sim.sysNL.plot_forcing(fig=fig, axes=axes, includeCK=False, label='Non linear', c='k')
    #sim.sysLI.plot_forcing(fig=fig, axes=axes, includeCK=True, plotCK0=True, label='Linear')
    #axes[0].legend()
plt.show()
