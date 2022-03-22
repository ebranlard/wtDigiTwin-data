"""

IMPORTANT:
   Remember that HydroStatic will give a moment (Mean vertical buoyancy force but also moment due to diameter)
   The moment will also affect the frequency.
   It can be set to more or less counteract the moment from gravity, but it's better to just put the HD structure out of the water!

"""

import matplotlib.pyplot as plt
import os
from numpy.linalg import inv

from welib.yams.models.simulator import *
from welib.yams.models.generator_oneRigidBody import generateOneRigidBodyModel
from welib.yams.models.generator import generateModel
np.set_printoptions(linewidth=300, precision=5)

from welib.FEM.utils import rigidTransformationTwoPoints, rigidTransformationTwoPoints_Loads
from welib.moor.mappp import Map




# --- Script parameters

figExport=True
create=True
# create=False
tMax = 10
tMax = None


# Model parameters
qop=None # <<<< VERY IMPORTANT FOR LIN MODEL
qdop=None # <<<< VERY IMPORTANT FOR LIN MODEL
CG_on_z =False # <<< Very important for most models with RNA

g=None
Refz=0
Moor=True

# --- Spar
fstFilename='Spar_MD1HD0SD0_F5T1N0S1/Main.fst'; modelName='F5T1N0S1_fnd'; z_BM=0; Refz=20; Moor=True; qop=[0,0.00,0,-0.0*np.pi/180, 0.0*np.pi/180,0,0]; aero_forces=False;  qdop=[0,0,0,0,0,0,10/60*2*np.pi];

# fstFilename='Spar_MD1HD0SD0_F100010/Main.fst'; modelName='B100010_moorO'; z_BM=0; Refz=20; Moor=True
# fstFilename='SparNoRNA_MD1F000010/Main.fst'; modelName='B000010_moorO'; z_BM=0; Refz=25; Moor=True
# fstFilename='SparNoRNA_MD0F000010/Main.fst'; modelName='B000010'; z_BM=0; Refz=25; Moor=True


# --- Tetra Spar
# fstFilename='Decays/Surge_MAP.fst'; modelName='B100000_moorM'; z_BM=0; Refz=16; Moor=True
# fstFilename='Decays/Pitch_MAP.fst'; modelName='B000010_moorM'; z_BM=0; Refz=16; Moor=True
# fstFilename='Decays/B100010_MAP.fst'; modelName='B100010_moorM'; z_BM=0; Refz=16; Moor=True; qop=[-0.08, -0.1720*np.pi/180]
# fstFilename='TS_MD1HD0SD0_F100010/Main.fst'; modelName='B100010_moorM'; z_BM=0; Refz=16; Moor=True; qop=[-0.08, -0.1720*np.pi/180]
# fstFilename='TS_MD1HD0SD0_F110111/Main.fst'; modelName='B110111_moorM'; z_BM=0; Refz=16; Moor=True; qop=[-0.08,0,0,-0.1720*np.pi/180, 0*np.pi/180]
# fstFilename='TS_MD1HD0SD0_F100010T0RNA/Main.fst'; modelName='F100010T0RNA_fnd'; z_BM=0; Refz=16; Moor=True; qop=[-0.08,-0.1720*np.pi/180]; aero_forces=False; 
# fstFilename='TS_MD1HD0SD0_F100010T1N0S1/Main.fst'; modelName='F100010T1N0S1_fnd'; z_BM=0; Refz=16; Moor=True; qop=[-0.0916,-0.0032, -0.1457,0]; aero_forces=False;  
# fstFilename='TS_MD1HD0SD0_F5T1N0S1/Main.fst'; modelName='F5T1N0S1_fnd'; z_BM=0; Refz=16; Moor=True; qop=[-0.076,0.00,0,-0.1830*np.pi/180, 0.042*np.pi/180,-0.1456,0]; aero_forces=False;  qdop=[0,0,0,0,0,0,10/60*2*np.pi];

# --- Tetra Spar Simple Map
# fstFilename='SimpleMAP/TetraSpar_TS_MD0F000010_SimpleMAP/Main.fst'; modelName='B000010'; z_BM=0 ; Moor=False
# fstFilename='SimpleMAP/TetraSpar_TS_MD1F100010_SimpleMAP/Main.fst'; modelName='B100010_moorM'; z_BM=0; Refz=25;

# --- No Ref
# fstFilename='SimpleMAP/SparNoRNA_TS_MD1F100000_NoRef_At0/Main.fst'; modelName='B100000_moorM'; z_BM=0;
# fstFilename='SimpleMAP/SparNoRNA_TS_MD1F100000_NoRef/Main.fst'; modelName='B100000_moorM'; z_BM=0;  Refz=0;
# fstFilename='SimpleMAP/SparNoRNA_TS_MD1F000010_NoRef/Main.fst'; modelName='B000010_moorM'; z_BM=0 ; Refz=0;
# fstFilename='SimpleMAP/SparNoRNA_TS_MD0F000010_NoRef/Main.fst'; modelName='B000010'; z_BM=0 ; Moor=False

# --- Simpler MAP
# fstFilename='SimpleMAP/SparNoRNA_TS_MD1F000010_SimpleMAP/Main.fst'; modelName='B000010_moorM'; z_BM=0; Refz=25;
# fstFilename='SimpleMAP/SparNoRNA_TS_MD0F000010_SimpleMAP/Main.fst'; modelName='B000010'; z_BM=0 ; Moor=False
# fstFilename='SimpleMAP/SparNoRNA_TS_MD1F100010_SimpleMAP/Main.fst'; modelName='B100010_moorM'; z_BM=0; Refz=25;

# --- Smaller mass
# fstFilename='SimpleMAP/SparNoRNA_TS_MD1F000010_NoRef_NewMass/Main.fst'; modelName='B000010_moorM'; z_BM=-10; g=0
# fstFilename='SimpleMAP/SparNoRNA_TS_MD0F000010_NoRef_NewMass/Main.fst'; modelName='B000010'; z_BM=0 ; Moor=False; g=10 
# fstFilename='SimpleMAP/SparNoRNA_TS_MD1F000010_NoRef_NewMass_NoHydroStat/Main.fst'; modelName='B000010_moorM'; z_BM=-10; g=10


# figSize =(10,4.2)
figSize =(12,10)
figTitle='StructMoorSim - ' + os.path.dirname(fstFilename)+' ' + modelName


P_HDRef = np.array((0,0,0))
P_EDRef = np.array((0,0,Refz))

T_ED2HD  = rigidTransformationTwoPoints(P_EDRef, P_HDRef)
T_HD2ED_l= rigidTransformationTwoPoints_Loads(P_HDRef, P_EDRef)


# --- MAP linearization
if Moor is True:
    moor = Map(fstFilename)
    # --- Plot initial 
#     fig, ax = moor.plot(numPoints=20, colors=['k'], ls='--')
    # --- Linearization with no displacement (NOTE: is it in equilibrium?)
    epsilon = 1e-2 # finite difference epsilon
    K_Moor = moor.linear(epsilon)    
    K_Moor = np.around(K_Moor,1)
else:
    K_Moor=np.zeros((6,6))

print("Mooring stiffness matrix (0,0,0)")
print(K_Moor)

K_Moor_ED1 = T_HD2ED_l.dot(K_Moor)
K_Moor_ED2 = T_HD2ED_l.dot(K_Moor.dot(T_ED2HD))
print('>>>T \n',T_ED2HD)
print('>>>TT\n',T_HD2ED_l)
print('>>>TT\n',T_HD2ED_l.T)

raise Exception()

print("Mooring stiffness matrix (0,0,Refz)")
print(K_Moor_ED1)

print("Mooring stiffness matrix (0,0,Refz)")
print(K_Moor_ED2)

K_Moor =K_Moor_ED2


k = K_Moor[0,0]

# --- Create symbolic model
if create:
    if modelName[0]=='B':
        generateOneRigidBodyModel(modelName, CG_on_z=CG_on_z)
    else:
        generateModel(modelName, aero_forces=False)

print('----------------------- SETUP SIMULATION -----------------------------------------')
WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)  # TODO
sim=SimulatorFromOF(WT, modelName=modelName, packageDir='py')
if modelName[0]=='B':
    time, dfFS, p = sim.setupSim(tMax=tMax, flavor='onebody', J_at_Origin=True)
else:
    time, dfFS, p = sim.setupSim(tMax=tMax, J_at_Origin=True)

sim.qop=qop
sim.qdop=qdop
if g is not None:
    sim.p['g'] = g

for i in range(6):
    for j in range(6):
        if j>=i:
            sim.p['KM_{}{}'.format(i,j)] = K_Moor[i,j]

sim.p['z_BM'] = z_BM
sim.p['z_T0'] = -Refz

sim.simulate(out=True, prefix='')
fig = sim.plot(export=figExport, figSize=figSize, title=figTitle)
p=sim.p




# --------------------------------------------------------------------------------}
# --- DEBUG OUT 
# --------------------------------------------------------------------------------{
# print('z_BG',p['z_BG'])
# print('x_BG',p['x_BG'])

# print('M    :',sim.M0)
# print('M_lin:',sim.M_lin)
# print('K_lin:',sim.K_lin)
# # print(sim.C_lin)
# # print(sim.B_lin)
# # 
# # print(sim.K_lin)
# J = p['J_yy_B']
# M = p['M_B']
# z_BG = p['z_BG']
# g    = p['g']
# print('J    :', J)
# print('M    :', M)
# print('zG   :', z_BG)
# print('g    :', g)
# Kg = -M*g*z_BG
# Km = z_BG**2 * k
# print('Kg   :',Kg)
# print('Km   :',Km)
# 
# om_g = np.sqrt(Kg/J)
# f_g  = om_g/(2*np.pi)
# T_g  = 1/f_g
# print('Gravity Only: Frequency = {:.2f} - T={:.1f}'.format(f_g, T_g))
# 
# om = np.sqrt((Kg+Km)/J)
# f  = om/(2*np.pi)
# T  = 1/f
# print('G+Mooring   : Frequency = {:.2f} - T={:.2f} (should be 11.04)'.format(f, T))

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
