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
from welib.fast.hydrodyn import HydroDyn
from welib.FEM.utils import rigidTransformationTwoPoints, rigidTransformationTwoPoints_Loads
np.set_printoptions(linewidth=300, precision=5)


modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010/Main.fst'; 
modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010_NoRef/Main.fst'; 
modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010_NoRef_NoCdCpCa/Main.fst'
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010_NoCdCpCa/Main.fst'; hydro=True 
fstFilenames=[
'SparNoRNA_MD0HD1SD0_F000010_NoRef_NoCdCpCa/Main.fst',
'SparNoRNA_MD0HD1SD0_F000010_NoCdCpCa/Main.fst']

fstFilenames=[
'TS_MD0HD1SD0_F101010_NoCdCpCa/Main.fst']


for  fstFilename in fstFilenames:
    print('----------------------- SETUP SIMULATION -----------------------------------------')
    print('>>> fstFileName: ',fstFilename)
    WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)  # TODO
    sim = SimulatorFromOF(WT, modelName=modelName, packageDir='py')
    time, dfFS, p = sim.setupSim(tMax=0, flavor='onebody', J_at_Origin=True)
    zRef = -sim.p['z_B0']
    print('>>>>>>>> zRefPtfm',zRef)

    # --- Linear Hydro
    q0=np.zeros(6)
    #     q0[4]=0.1*np.pi/180
    hd = HydroDyn(fstFilename)

    # ---
    MCKFh = hd.linearize_RigidMotion2Loads(q0, RefPointMotion=(0,0,zRef), RefPointMapping=(0,0,zRef))
    Kh = pd.DataFrame(data=MCKFh[2], columns=['x','y','z','phi_x','phi_y','phi_z'], index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])
    print('>>> Kh\n',np.around(Kh,1))
    # ---
    MCKFh = hd.linearize_RigidMotion2Loads(q0, RefPointMotion=(0,0,0), RefPointMapping=(0,0,zRef))
    Kh = pd.DataFrame(data=MCKFh[2], columns=['x','y','z','phi_x','phi_y','phi_z'], index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])
    print('>>> Kh\n',np.around(Kh,1))
    # ---
    MCKFh = hd.linearize_RigidMotion2Loads(q0, RefPointMotion=(0,0,0), RefPointMapping=(0,0,0))
    Kh = pd.DataFrame(data=MCKFh[2], columns=['x','y','z','phi_x','phi_y','phi_z'], index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])
    print('>>> Kh\n',np.around(Kh,1))
    # ---
    MCKFh = hd.linearize_RigidMotion2Loads(q0, RefPointMotion=(0,0,zRef), RefPointMapping=(0,0,0))
    Kh = pd.DataFrame(data=MCKFh[2], columns=['x','y','z','phi_x','phi_y','phi_z'], index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])
    print('>>> Kh\n',np.around(Kh,1))

    su = sim.pkg.info()['su']
    sq = sim.WT.DOFname

#     if zRef==0:
#         zRef = 25
#         P_HDRef   = np.array((0,0,0))
#         P_EDRef   = np.array((0,0,zRef))
#         T_HD2ED   = rigidTransformationTwoPoints(P_HDRef, P_EDRef)
#         T_ED2HD   = rigidTransformationTwoPoints(P_EDRef, P_HDRef)
#         T_HD2ED_l = rigidTransformationTwoPoints_Loads(P_HDRef, P_EDRef)
#         T_ED2HD_l = rigidTransformationTwoPoints_Loads(P_EDRef, P_HDRef)
#         print(T_HD2ED_l)
#         print(T_HD2ED)
#         
#         Kh2 = T_HD2ED_l.dot(Kh.dot(T_HD2ED)) # Transform from Refz to 0
#         Kh2 = pd.DataFrame(data = Kh2, columns = Kh.columns, index=Kh.index)
#         print('Kh2\n',Kh2)

# 
# 
# Kh2 = T_ED2HD_l.dot(Kh) # Transform from Refz to 0
# Kh2 = pd.DataFrame(data = Kh2, columns = Kh.columns, index=Kh.index)
# 
# Kh_ = Kh.loc[su,sq]
# print('>>> Kh_\n',Kh_)
# 
# Kh2_ = Kh2.loc[su,sq]
# print('>>> Kh2_\n',Kh2_)
