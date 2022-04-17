""" 

NOTES:
- Results are fine with no drag and ref=0
    - TODO problem with ref
    - TODO deal with Drag

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

# ---- Script parameters
create=False
# create=True
runSim=True

tMax = 5
tMax = 20
tMax = None

Method = 1
qop=None
qdop=None

# --- Spar No RNA
# modelName = 'B001000_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F001000/Main.fst'; 
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010_NoRef/Main.fst';  
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010_NoRef_NoCdCpCa/Main.fst';  
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010_NoCd/Main.fst';  
# modelName = 'B000010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F000010/Main.fst';  
# modelName = 'B101010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F101010_NoRef/Main.fst'; 
# modelName = 'B101010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F101010/Main.fst'; 
# modelName = 'B111111_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F111111_NoRef/Main.fst'; 
# modelName = 'B111111_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F111111_NoCd/Main.fst'; 
# modelName = 'B111111_hydroO'; fstFilename = 'SparNoRNA_MD0HD1SD0_F111111_NoCd/Main.fst'; 
# modelName = 'B111111_hydroO'; fstFilename = 'SparNoRNA_MD0HD1SD0_F111111/Main.fst'; 

# --- Spar
modelName = 'B001000_hydro0'; fstFilename = 'Spar_MD0HD1SD0_F001000/Main.fst';  qop=[1.20] # OK
# modelName = 'B111111_hydro0'; fstFilename = 'Spar_MD0HD1SD0_F111111_NoCd/Main.fst';  qop=[-1.5, 0.8, 1.23, 0,0,0.14*np.pi/180]
# modelName = 'B111111_hydroO'; fstFilename = 'Spar_MD0HD1SD0_F111111_NoCd/Main.fst';  qop=[-1.6, 0.8, 1.20, 0,0.0*np.pi/180,0.14*np.pi/180]# NOTE: qop Pitch will introduce instability
# modelName = 'B111111_hydro0'; fstFilename = 'Spar_MD0HD1SD0_F111111/Main.fst';       qop=[-1.6, 0.8, 1.23, 0,0,0.14*np.pi/180]
# modelName = 'B111111_hydroO'; fstFilename = 'Spar_MD0HD1SD0_F111111/Main.fst';       qop=[-1.6, 0.8, 1.23, 0,0,0.14*np.pi/180]

# --- Tetra Spar
# modelName = 'B000010_hydro0'; fstFilename = 'TS_MD0HD1SD0_F000010_NoCdCpCa/Main.fst'; qop=[-2.266*np.pi/180];
# modelName = 'B101010_hydro0'; fstFilename = 'TS_MD0HD1SD0_F101010/Main.fst'; 
# modelName = 'B101010_hydro0'; fstFilename = 'TS_MD0HD1SD0_F101010_NoCdCpCa/Main.fst';  qop=[-2.5, 3.2, -1.68*np.pi/180]
# modelName = 'B101010_hydro0'; fstFilename = 'TS_MD0HD1SD0_F101010_NoCd/Main.fst';  qop=[-2.5, 3.2, -1.68*np.pi/180]
# modelName = 'B101010_hydroO'; fstFilename = 'TS_MD0HD1SD0_F101010_NoCd/Main.fst';  qop=[-2.5, 3.2, -1.68*np.pi/180]
# modelName = 'B111111_hydro0'; fstFilename = 'TS_MD0HD1SD0_F111111_NoCdCpCa/Main.fst'; 
# modelName = 'B111111_hydro0'; fstFilename = 'TS_MD0HD1SD0_F111111_NoCd/Main.fst';   qop=[-2.8, 0.08, 3.2, 0, -2.0*np.pi/180, 0.000*np.pi/180]
# modelName = 'B111111_hydroO'; fstFilename = 'TS_MD0HD1SD0_F111111_NoCd/Main.fst';   qop=[-2.8, 0.08, 3.2, 0, -2.0*np.pi/180, 0.000*np.pi/180]


# --- Generate python package
if create:
    if modelName[0]=='B':
        generateOneRigidBodyModel(modelName)
    else:
        generateModel(modelName, aero_forces=False, moor_loads=False, hydro_loads=True)

MCKh=None




# MCKh=0
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
    zRef = -sim.p['z_B0']
    su = sim.pkg.info()['su']
    sq = sim.WT.DOFname

    # --- Linear Hydro
    print('----------------------- LINEAR HYDRO  --------------------------------------------')
    q0=np.zeros(6)
    if qop is not None:
        q0_ = pd.DataFrame(data=q0, index=['x','y','z','phi_x','phi_y','phi_z'])
        for i,s in enumerate(sq): 
            q0_.loc[s] =qop[i]
        q0 = q0_.values.flatten()
    hd = HydroDyn(fstFilename)
    if MCKh is 0:
        Mh=np.zeros((6,6))
        Ch=np.zeros((6,6))
        Kh=np.zeros((6,6))
    if MCKh is None:
        if 'hydroO' in modelName:
            MCKFh = hd.linearize_RigidMotion2Loads(q0, RefPointMotion=(0,0,zRef), RefPointMapping=(0,0,zRef) ) # <<< Good if hydroO model
        else:
            MCKFh = hd.linearize_RigidMotion2Loads(q0, RefPointMotion=(0,0,zRef), RefPointMapping=(0,0,0) ) # <<< Good if hydro0 model
#       MCKFh = hd.linearize_RigidMotion2Loads(q0, RefPointMotion=(0,0,0), RefPointMapping=(0,0,0) ) # OLD and BAD
        Mh,Ch,Kh,Fh0=MCKFh
    print('Ch\n',Ch)

    # TODO this below only works for subset
    Mh = pd.DataFrame(data=Mh, columns=['x','y','z','phi_x','phi_y','phi_z'], index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])
    Ch = pd.DataFrame(data=Ch, columns=['x','y','z','phi_x','phi_y','phi_z'], index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])
    Kh = pd.DataFrame(data=Kh, columns=['x','y','z','phi_x','phi_y','phi_z'], index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])
    Fh = pd.DataFrame(data=Fh0, index=['F_hx','F_hy','F_hz','M_hx','M_hy','M_hz'])
#     Refz = - sim.p['z_B0']
#     P_HDRef   = np.array((0,0,0))
#     P_EDRef   = np.array((0,0,Refz))
#     T_HD2ED   = rigidTransformationTwoPoints(P_HDRef, P_EDRef)
#     T_ED2HD   = rigidTransformationTwoPoints(P_EDRef, P_HDRef)
#     T_HD2ED_l = rigidTransformationTwoPoints_Loads(P_HDRef, P_EDRef)
#     T_ED2HD_l = rigidTransformationTwoPoints_Loads(P_EDRef, P_HDRef)
# 
#     # Mh_ = T_HD2ED_l.dot(K_Moor.dot(T_ED2HD)) # Transform from 0 to Refz
#     Mh2 = T_ED2HD_l.dot(Mh.dot(T_HD2ED)) # Transform from Refz to 0
#     Kh2 = T_ED2HD_l.dot(Kh.dot(T_HD2ED)) # Transform from Refz to 0
#     Ch2 = T_ED2HD_l.dot(Ch.dot(T_HD2ED)) # Transform from Refz to 0
#     Mh2 = pd.DataFrame(data = Mh2, columns = Mh.columns, index=Mh.index)
#     Kh2 = pd.DataFrame(data = Kh2, columns = Kh.columns, index=Kh.index)
#     Ch2 = pd.DataFrame(data = Ch2, columns = Ch.columns, index=Ch.index)
# 
# 
#     Mh2 = T_ED2HD_l.dot(Mh) # Transform from Refz to 0
#     Kh2 = T_ED2HD_l.dot(Kh) # Transform from Refz to 0
#     Ch2 = T_ED2HD_l.dot(Ch) # Transform from Refz to 0
#     Mh2 = pd.DataFrame(data = Mh2, columns = Mh.columns, index=Mh.index)
#     Kh2 = pd.DataFrame(data = Kh2, columns = Kh.columns, index=Kh.index)
#     Ch2 = pd.DataFrame(data = Ch2, columns = Ch.columns, index=Ch.index)
# 
#     print('>>>>>>>>>> HACK DAMPING <<<<<')
#     Ch.loc['F_hy','y']    =100000000
#     Ch.loc['M_hz','phi_z']=100000000
#     print('>>>>>>>>>> HACK STIFFNESS <<<<<')
#     Kh.loc['F_hy','y']    +=1000000
#     Kh.loc['M_hz','phi_z']+=1000000

    Mh_ = Mh.loc[su,sq]
    Ch_ = Ch.loc[su,sq]
    Kh_ = Kh.loc[su,sq]
    Fh_ = Fh.loc[su]
    print('>>> Ch_\n',Ch_)
    print('>>> Mh_\n',Mh_)
    print('>>> Kh_\n',Kh_)
    print('>>> Fh_\n',Fh_)
# 
#     Mh2_ = Mh2.loc[su,sq]
#     Ch2_ = Ch2.loc[su,sq]
#     Kh2_ = Kh2.loc[su,sq]
#     print('>>> Ch2_\n',Ch2_)
#     print('>>> Mh2_\n',Mh2_)
#     print('>>> Kh2_\n',Kh2_)
# #     Ch_ = Ch_*0.2
# #     Mh_ = Mh_*0
# #     Kh_ = Kh_*0
    MCKu = Mh_, Ch_, Kh_
#     MCKu = Mh2_, Ch2_, Kh2_

#     print('>>>>>>>>> zB0',sim.p['z_B0'])
#     sim.p['z_B0']=0

    # --- For non linear simulation (NOTE: need special solver with qdd...)
    #     umesh = hd.u['Morison']['Mesh']
    #     ymesh = hd.y['Morison']['Mesh']
    #             # Reference point motion
    #             q  = dfOF[qCol].iloc[it].values
    #             qd = dfOF[qdCol].iloc[it].values
    #             omega = (qd[3],qd[4],qd[5]) # ...
    #             # Rigid body motion of the mesh
    #             umesh.rigidBodyMotion(q=q, qd=qd, qdd=qdd)
    #             # Calculate hydrodynamic loads at every nodes 
    #             hd.y=hd.calcOutput(t, u=hd.u, y=hd.y, optsM=optsM)
    #             # Store mesh
    #             msy.store(ymesh, it)
    #             # Compute integral loads (force&moment) at the reference point (translated but not rotated)
    #             fh,Mh = ymesh.mapLoadsToPoint((q[0],q[1],q[2]))
    #     #sim.setCoupledHydroLoads()
    #     if Method==1:
    # 
    #     elif Method==2:
    #         MCKu=None
    #         MM0      = pkg.mass_matrix(q0,p)
    #         forcing0 = pkg.forcing(t,q0,qd0,p,u)
    #         # --- integrate non-linear system
    #         fM = lambda x: pkg.mass_matrix(x, p)
    #         fF = lambda t,x,xd: pkg.forcing(t, x, xd, p=p, u=u)
    #         sysNL = MechSystem(fM, F=fF, x0=q0, xdot0=qd0 )

    # --- uop
    print('----------------------- OPERATING POINT ------------------------------------------')
    print('>>>> UOP')
    uop = sim.uop
#     for s in su:
#         uop[s]=Fh_.loc[s].values[0]
#     uop['F_hz'] = p['M_B']*p['g']
    print(uop['F_hz'])
    print(p['M_B']*p['g'])
    print(uop)

    #uop=None
    #du=None
    sim.qop  = qop
    sim.qdop = qdop
    # Using mean as op
    # qop  = np.array([np.mean(dfFS[c]) for c in WT.q_channels])
    # qdop = np.array([np.mean(dfFS[c]) for c in WT.qd_channels])*0
    # --- Simulation
    #sim.setInputs(u, du, uop, qop, qdop)
    sim.simulate(out=True, prefix='_hydroPyCoupled_'+modelName, MCKu=MCKu, NL=False)
    sim.plot(export=True,  prefix='_hydroPyCoupled_'+modelName, nPlotCols=3)

    # --- Plot forcing
    #fig,axes = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    #fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    #axes=np.atleast_1d(axes)
    #sim.sysNL.plot_forcing(fig=fig, axes=axes, includeCK=False, label='Non linear', c='k')
    #sim.sysLI.plot_forcing(fig=fig, axes=axes, includeCK=True, plotCK0=True, label='Linear')
    #axes[0].legend()
plt.show()
