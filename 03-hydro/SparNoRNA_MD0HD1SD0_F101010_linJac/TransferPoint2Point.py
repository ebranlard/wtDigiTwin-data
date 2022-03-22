""" 
Convert from RefPoint (0,0,0) to Ptfm and back

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from welib.FEM.utils import rigidTransformationTwoPoints
from welib.yams.rotations import BodyXYZ_A, SmallRot_DCM

Refz=16

P1 = np.array((0,0,0))
P2 = np.array((0,0,Refz))

# Rigid transformation matrix between DOFs of node j and k where node j is the leader node (Pk-Pj)
T= rigidTransformationTwoPoints(P2, P1)
print(T)



qCol1  = ['PRPSurge_[m]'    ,'PRPSway_[m]'    ,'PRPHeave_[m]'   ,'PRPRoll_[rad]'    ,'PRPPitch_[rad]'   ,'PRPYaw_[rad]'     ]
qdCol1 = ['PRPTVxi_[m/s]'   ,'PRPTVyi_[m/s]'  ,'PRPTVzi_[m/s]'  ,'PRPRVxi_[rad/s]'  ,'PRPRVyi_[rad/s]'  ,'PRPRVzi_[rad/s]'  ]
qddCol1= [ 'PRPTAxi_[m/s^2]','PRPTAyi_[m/s^2]','PRPTAzi_[m/s^2]','PRPRAxi_[rad/s^2]','PRPRAyi_[rad/s^2]','PRPRAzi_[rad/s^2]']
qCol2  = ['Q_Sg_[m]'   ,'Q_Sw_[m]'   ,'Q_Hv_[m]'   ,'Q_R_[rad]'   ,'Q_P_[rad]'   ,'Q_Y_[rad]']
qdCol2 = ['QD_Sg_[m/s]','QD_Sw_[m/s]','QD_Hv_[m/s]','QD_R_[rad/s]','QD_P_[rad/s]','QD_Y_[rad/s]']
qddCol2= [None]*6




# --- Read OpenFAST
fstFilename = '../models/MD0HD1SD0_PitchDecayRigid_Heave/Main.fst'
dfOF = weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
time=dfOF['Time_[s]'].values

# --- Convert motion from Ref to Ptfm
#fh = np.zeros((len(time), 6))
qPtfm = np.zeros((len(time), 6))
qRef = np.zeros((len(time), 6))
qdPtfm = np.zeros((len(time), 6))
qdRef = np.zeros((len(time), 6))
qddPtfm = np.zeros((len(time), 6))
qddRef = np.zeros((len(time), 6))
for it, t in enumerate(time):
    # Reference point motion
    q   = dfOF[qCol1].iloc[it].values
    qd  = dfOF[qdCol1].iloc[it].values
    qdd = dfOF[qddCol1].iloc[it].values

    # Transform to platform
    r_AB0           = P2-P1
    omega           = qd[3:]
    omega_dot       = qdd[3:]
    R_b2g           = SmallRot_DCM(q[3], q[4], q[5]).T
    r_AB            = R_b2g.dot(r_AB0)
    om_x_r          = (np.cross(omega, r_AB))
    qPtfm  [it, :3] = q  [:3] + (r_AB - r_AB0)
    qdPtfm [it, :3] = qd [:3] + om_x_r
    qddPtfm[it, :3] = qdd[:3] + np.cross(omega_dot, r_AB) + np.cross(omega, om_x_r)
    qPtfm  [it,3:6] = q[3:6]
    qdPtfm [it,3:6] = omega
    qddPtfm[it,3:6] = omega_dot

    # Platform point motion
    q   = dfOF[qCol2].iloc[it].values
    qd  = dfOF[qdCol2].iloc[it].values
    #qdd = dfOF[qddCol2].iloc[it].values
    qdd = qddPtfm[it,:]

    # Transform to Ref
    r_AB0           = P1-P2
    omega           = qd[3:]
    omega_dot       = qdd[3:]
    R_b2g           = SmallRot_DCM(q[3], q[4], q[5]).T
    r_AB            = R_b2g.dot(r_AB0)
    om_x_r          = (np.cross(omega, r_AB))
    qRef   [it, :3] = q  [:3] + (r_AB - r_AB0)
    qdRef  [it, :3] = qd [:3] + om_x_r
    qddRef [it, :3] = qdd[:3] + np.cross(omega_dot, r_AB) + np.cross(omega, om_x_r)
    qRef   [it,3:6] = q[3:6]
    qdRef  [it,3:6] = omega
    qddRef [it,3:6] = omega_dot


    # Forces
    #fh[it,:] = -M.dot(qdd) - C.dot(qd) - K.dot(q)


# lCols = ['HydroFxi_[N]','HydroFyi_[N]','HydroFzi_[N]','HydroMxi_[N-m]','HydroMyi_[N-m]','HydroMzi_[N-m]']
# dfPH = pd.DataFrame(data=np.column_stack((time,fh)), columns=['Time_[s]']+lCols)


# --- Plot
fig,axes = plt.subplots(6, 3, sharey=False, figsize=(12.8,8.5)) # (6.4,4.8)
fig.subplots_adjust(left=0.11, right=0.95, top=0.95, bottom=0.07, hspace=0.40, wspace=0.22)
# DOF
for iCol, (col1,col2) in enumerate(zip(qCol1, qCol2)):
    # Positions
    axes[iCol,0].plot(time, dfOF[col1].values  , '-', label='Hydro Ref')
    axes[iCol,0].plot(time, dfOF[col2].values  , '--', label='Ptfm')
    axes[iCol,0].plot(time, qPtfm[:,iCol]      , ':',   label='Ptfm (Transfered)')
    axes[iCol,0].plot(time, qRef [:,iCol]      , '-.',  label='Hydro ref (Transfered)')
    axes[iCol,0].set_ylabel(col1.replace('_',' '))
for iCol, (col1,col2) in enumerate(zip(qdCol1, qdCol2)):
    # Velocities
    axes[iCol,1].plot(time, dfOF[col1].values  , '-', label='Hydro Ref')
    axes[iCol,1].plot(time, dfOF[col2].values  , '--', label='Ptfm')
    axes[iCol,1].plot(time, qdPtfm[:,iCol]      , ':',   label='Ptfm (Transfered)')
    axes[iCol,1].plot(time, qdRef [:,iCol]      , '-.',  label='Hydro ref (Transfered)')
    axes[iCol,1].set_ylabel(col1.replace('_',' '))

for iCol, (col1,col2) in enumerate(zip(qddCol1, qddCol2)):
    # Acceleration
    axes[iCol,2].plot(time, dfOF[col1].values  , '-', label='Hydro Ref')
#     axes[iCol,2].plot(time, dfOF[col2].values  , '--', label='Ptfm')
    axes[iCol,2].plot(time, qddPtfm[:,iCol]      , ':',   label='Ptfm (Transfered)')
    axes[iCol,2].plot(time, qddRef [:,iCol]      , '-.',  label='Hydro ref (Transfered)')
    axes[iCol,2].set_ylabel(col1.replace('_',' '))


# # Forces
# for iCol, col in enumerate(lCols):
#     axes[iCol,1].plot(time, dfPH[col].values+F0[iCol]  , label='Python linear')
#     axes[iCol,1].plot(time, dfOF[col].values           , 'k:', label='OpenFAST')
#     axes[iCol,1].set_ylabel(col.replace('_',' '))
axes[0,0].legend()
axes[5,0].set_xlabel('Time [s]')
axes[5,1].set_xlabel('Time [s]')

plt.show()

if __name__ == '__main__':
    pass
