import matplotlib.pyplot as plt
import os

from welib.yams.models.simulator import *
from welib.yams.models.generator_oneRigidBody import generateOneRigidBodyModel
from welib.yams.models.generator import generateModel
np.set_printoptions(linewidth=300, precision=5)

# --- Script parameters

figExport=True
# figExport=False
create=True
create=False
tMax = 10
tMax = None


# Model parameters
qop=None # <<<< VERY IMPORTANT FOR LIN MODEL
qdop=None # <<<< VERY IMPORTANT FOR LIN MODEL
CG_on_z =False # <<< Very important for most models with RNA

#--- Spar
fstFilename = 'Spar_F5T1N0S1/Main_Spar_ED.fst'; modelName='F5T1N0S1_fnd'; qop=[0,0,0,0,0,0,0]; qdop=[0,0,0,0,0,0,10/60*2*np.pi];

# --- Tetra Spar
# 
fstFilename = 'MD0HD0SD0_F5T1N0S1/Main.fst'   ; modelName='F5T1N0S1_fnd'; qop=[0.25,-0.5162, 0, -0.1865*np.pi/180, 0.27*np.pi/180 ,-0.14,0]; qdop=[0,0,0,0,0,0,10/60*2*np.pi]; # NOTE: pitch init is very important. Yaw varies a lot in full sim, so won't be able to capture
# fstFilename = 'MD0HD0SD0_F110111/Main.fst'   ; modelName='B110111'; qop=[0.34,-0.51, 0, -0.185*np.pi/180, 0.02*np.pi/180]

# --- Tetra Spar simple models
# fstFilename = 'MD0HD0SD0_F100010/Main.fst'   ; modelName='B100010'; #qop=[0,-0.155*np.pi/180]
# fstFilename = 'MD0HD0SD0_F100010/Main.fst'   ; modelName='F100010T0RNA_fnd_noLoads'; qop=[0,-0.155*np.pi/180]

# fstFilename = 'MD0HD0SD0_F000010/Main.fst'   ; modelName='B000010';  qop=[-0.155*np.pi/180] # TODO Auto
# fstFilename = 'MD0HD0SD0_F000010/Main.fst'   ; modelName='F000010T0RNA_fnd_noLoads'; qop=[-0.155*np.pi/180]


figTitle='StructSim - ' + os.path.dirname(fstFilename)+' ' + modelName

if create:
    if modelName[0]=='B':
        generateOneRigidBodyModel(modelName, CG_on_z=CG_on_z)
    else:
        generateModel(modelName, aero_forces=False, moor_loads=False, hydro_loads=False)

print('----------------------- SETUP SIMULATION -----------------------------------------')
WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)  # TODO
sim=SimulatorFromOF(WT, modelName=modelName, packageDir='py')
if modelName[0]=='B':
    time, dfFS, p = sim.setupSim(tMax=tMax, flavor='onebody', J_at_Origin=True)
else:
    time, dfFS, p = sim.setupSim(tMax=tMax, J_at_Origin=True)

sim.qop=qop
sim.qdop=qdop
sim.simulate(out=True, prefix='')
fig = sim.plot(export=figExport, figSize=(12,10), title=figTitle)
p=sim.p

# --------------------------------------------------------------------------------}
# --- DEBUG 
# --------------------------------------------------------------------------------{
# print('z_BG',p['z_BG'])
# print('x_BG',p['x_BG'])

# print(sim.M0)
# print(sim.M_lin)
# print(sim.K_lin)
# print(sim.C_lin)
# print(sim.B_lin)
# 
# print(sim.K_lin)
# 
# print('q0          ', sim.q0          )
# print('q0_lin      ', sim.q0_lin      )
# print('forcing0    ', sim.forcing0    )
# print('forcing0_lin', sim.forcing0_lin)
# print('qop         ', sim.qop         )
# print('qdop        ', sim.qdop        )
# print('uop         ', sim.uop         )
# F100010
# [[ 5.58417e+06 -3.06176e+08]
#  [-3.06176e+08  2.66767e+10]]
# [[ 5.58417e+06 -3.06240e+08]
#  [-3.06240e+08  2.66767e+10]]
# [[0.00000e+00 0.00000e+00]
#  [0.00000e+00 3.00319e+09]]
# [[0. 0.]
#  [0. 0.]]
# [[ 0.99144]
#  [69.58922]]
# q0           [0.      0.01745 0.      0.     ]
# q0_lin       [0.      0.01745 0.      0.     ]
# forcing0     [        0.      -62347987.35301]
# forcing0_lin [0. 0.]
# qop          None
# qdop         None
# uop          {'T_a': 0}

# --- B100010
# [[ 5.58417e+06 -3.06194e+08]
#  [-3.06194e+08  2.66767e+10]]
# [[ 5.58417e+06 -3.06240e+08]
#  [-3.06240e+08  2.66767e+10]]
# [[0.00000e+00 0.00000e+00]
#  [0.00000e+00 3.00319e+09]]
# [[0. 0.]
#  [0. 0.]]
# [[0. 0.]
#  [0. 0.]]
# q0           [0.      0.01745 0.      0.     ]
# q0_lin       [0.      0.01745 0.      0.     ]
# forcing0     [        0.      -52412899.09496]


def main():
    # --- Rigid "F2"
    #fstFilename = '_F2T0RNANoRefH/Main_Spar_ED.fst' ;modelName='F2T0RNA_fnd';sim_name='F2T0RNA_NoRefH'
    #fstFilename = '_F2T0_NoRNA_NoRefH/Main_Spar_ED.fst' ;modelName='F2T0RNA_fnd';sim_name='F2T0RNA_NoRNA_NoRefH'
    #fstFilename = '_F2T0RNA/Main_Spar_ED.fst'      ;modelName='F2T0RNA_fnd';sim_name='F2T0RNA'

    #fstFilename = '_F2T0N0S1/Main_Spar_ED.fst'; modelName='F2T0N0S1_fnd'; sim_name='F2T0N0S1'

    # Phi x phi z
    #fstFilename = '_F000101T0RNA/Main_Spar_ED.fst'; modelName='B000101'; sim_name=modelName
    #fstFilename = '_F000101T0RNA/Main_Spar_ED.fst'; modelName='F000101T0RNA_fnd'; sim_name=modelName
    #fstFilename = '_F000101T0N0S1/Main_Spar_ED.fst'; modelName='F000101T0N0S1_fnd'; sim_name=modelName

    # Phi y phi z
    #fstFilename = '_F000011T0RNA/Main_Spar_ED.fst'; modelName='B000011'; sim_name=modelName

    # Phi x phi y phi z 
    #fstFilename = '_F000111T0RNA/Main_Spar_ED.fst'; modelName='B000111'; sim_name=modelName
    #fstFilename = '_F000111T0RNA/Main_Spar_ED.fst'; modelName='F000111T0RNA_fnd'; sim_name=modelName

    # --- Rotor
    #fstFilename = '_F000111T0N0S1/Main_Spar_ED.fst'; modelName='F000111T0N0S1_fnd'; sim_name=modelName
    #fstFilename = '_F5T0N0S1/Main_Spar_ED.fst'; modelName='F5T0N0S1_fnd'; sim_name='F5T0N0S1'

    # --- Flexibility "T1, T2"
    #fstFilename = '_F0T1RNA/Main_Spar_ED.fst'; modelName='F0T1RNA'; sim_name='F0T1RNA'

    #fstFilename = '_F0T2_NoRNA_sym/Main_Spar_ED.fst'; modelName='F0T2RNA'; sim_name='F0T2_NoRNA_sym'  # NOTE: Works fine large disp, symmetric shapes, with HubMass and NacMass, Twr2Shaft, detoriate slightly with overhang 
    #fstFilename = '_F0T2_NoRNA/Main_Spar_ED.fst'; modelName='F0T2RNA'; sim_name='F0T2_NoRNA'  # NOTE: with asymmetric shape functions, cannot achieve as good a result somehow. Wrong alpha???

    #fstFilename = '_F0T2RNA/Main_Spar_ED.fst'; modelName='F0T2RNA'; sim_name='F0T2RNA'
    #fstFilename = '_F0T2RNA_sym/Main_Spar_ED.fst'; modelName='F0T2RNA'; sim_name='F0T2RNA_sym'

    #fstFilename = '_F0T2N0S1/Main_Spar_ED.fst'; modelName='_F0T2N0S1'; sim_name=modelName;

    # --- Floater + Flexibility "F2T1"
    #fstFilename = '_F2T1RNANoRefH/Main_Spar_ED.fst'; modelName='F2T1RNA_fnd'; sim_name='F2T1RNA_NoRefH'
    #fstFilename = '_F2T1RNA_SmallAngle/Main_Spar_ED.fst'; modelName='F2T1RNA_fnd'; sim_name='F2T1RNA_SmallAngle'
    #fstFilename = '_F2T1RNA/Main_Spar_ED.fst'; modelName='F2T1RNA_fnd'; sim_name='F2T1RNA_LargeAngle'

    # --- "F3 T1"
    #fstFilename = '_F3T1RNA/Main_Spar_ED.fst'; modelName='F3T1RNA_fnd'; sim_name='F3T1RNA'
    # --- "F5 T1"
    pass
    # --- Print parameters
#     print('--------------------')
#     print('Strucural Parameters: ')
#     for k,v in p.items():
#         if hasattr(v,'__len__'):
#             print('{:10s}:\n{}'.format(k,v))
#         else:
#             print('{:10s}:{}'.format(k,v))

if __name__ == '__main__':
    plt.show()
