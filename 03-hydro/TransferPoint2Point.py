""" 
Convert from RefPoint (0,0,0) to Ptfm and back

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio
from welib.FEM.utils import rigidTransformationTwoPoints, rigidTransformationTwoPoints_Loads
from welib.yams.kinematics import rigidBodyMotion2Points_q6

Refz=16

P_HDRef = np.array((0,0,0))
P_EDRef = np.array((0,0,Refz))




qCol1  = ['PRPSurge_[m]'    ,'PRPSway_[m]'    ,'PRPHeave_[m]'   ,'PRPRoll_[rad]'    ,'PRPPitch_[rad]'   ,'PRPYaw_[rad]'     ]
qdCol1 = ['PRPTVxi_[m/s]'   ,'PRPTVyi_[m/s]'  ,'PRPTVzi_[m/s]'  ,'PRPRVxi_[rad/s]'  ,'PRPRVyi_[rad/s]'  ,'PRPRVzi_[rad/s]'  ]
qddCol1= [ 'PRPTAxi_[m/s^2]','PRPTAyi_[m/s^2]','PRPTAzi_[m/s^2]','PRPRAxi_[rad/s^2]','PRPRAyi_[rad/s^2]','PRPRAzi_[rad/s^2]']
qCol2  = ['Q_Sg_[m]'     , 'Q_Sw_[m]'     , 'Q_Hv_[m]'     , 'Q_R_[rad]'     , 'Q_P_[rad]'     , 'Q_Y_[rad]']
qdCol2 = ['QD_Sg_[m/s]'  , 'QD_Sw_[m/s]'  , 'QD_Hv_[m/s]'  , 'QD_R_[rad/s]'  , 'QD_P_[rad/s]'  , 'QD_Y_[rad/s]']
qddCol2= ['QD2_Sg_[m/s^2]' , 'QD2_Sw_[m/s^2]' , 'QD2_Hv_[m/s^2]' , 'QD2_R_[rad/s^2]' , 'QD2_P_[rad/s^2]' , 'QD2_Y_[rad/s^2]']

lCols = ['HydroFxi_[N]','HydroFyi_[N]','HydroFzi_[N]','HydroMxi_[N-m]','HydroMyi_[N-m]','HydroMzi_[N-m]']



# --- Read OpenFAST
fstFilename = 'MD0HD1SD0_PitchDecayRigid_Heave/Main.fst'
dfOF = weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
time=dfOF['Time_[s]'].values


# NOTE: In theory these should depend on the operating point (q, qd, qdd)
T_HD2ED= rigidTransformationTwoPoints(P_HDRef, P_EDRef)
T_ED2HD= rigidTransformationTwoPoints(P_EDRef, P_HDRef)

T_HD2ED_l= rigidTransformationTwoPoints_Loads(P_HDRef, P_EDRef)
T_ED2HD_l= rigidTransformationTwoPoints_Loads(P_EDRef, P_HDRef)

# --- Convert motion from Ref to Ptfm
fh_HD = np.zeros((len(time), 6))
fh_ED = np.zeros((len(time), 6))
fh_HD2 = np.zeros((len(time), 6))
qRef        = np.zeros((len(time),6))
qdRef       = np.zeros((len(time),6))
qddRef      = np.zeros((len(time),6))
qPtfm       = np.zeros((len(time),6))
qdPtfm      = np.zeros((len(time),6))
qddPtfm     = np.zeros((len(time),6))
qPtfm_lin   = np.zeros((len(time),6))
qdPtfm_lin  = np.zeros((len(time),6))
qddPtfm_lin = np.zeros((len(time),6))
for it, t in enumerate(time):
    # HD Reference point motion
    q   = dfOF[qCol1].iloc[it].values
    qd  = dfOF[qdCol1].iloc[it].values
    qdd = dfOF[qddCol1].iloc[it].values
    # Transform to ED Ref platform
    qPtfm[it,:],qdPtfm[it,:], qddPtfm[it,:] = rigidBodyMotion2Points_q6(PSource0=P_HDRef, PDest0=P_EDRef, q=q, qd=qd, qdd=qdd)

    # Linear Rigid transformation matrix between DOFs of node j and k where node j is the leader node (Pk-Pj)
    qPtfm_lin[it,:]   = T_HD2ED.dot(q)
    qdPtfm_lin[it,:]  = T_HD2ED.dot(qd)
    qddPtfm_lin[it,:] = T_HD2ED.dot(qdd)


    # ED Platform point motion
    q   = dfOF[qCol2].iloc[it].values
    qd  = dfOF[qdCol2].iloc[it].values
    qdd = dfOF[qddCol2].iloc[it].values
    # Transform to HD Ref
    qRef[it,:],qdRef[it,:], qddRef[it,:] = rigidBodyMotion2Points_q6(PSource0=P_EDRef, PDest0=P_HDRef, q=q, qd=qd, qdd=qdd)

    # --- Loads
    fh_HD [it,:] = dfOF[lCols].iloc[it].values

    fh_ED [it,:] = T_HD2ED_l.dot(fh_HD[it,:])
    fh_HD2[it,:] = T_ED2HD_l.dot(fh_ED[it,:])
    #fh[it,:] = -M.dot(qdd) - C.dot(qd) - K.dot(q)
    #fh[it,:] = -M.dot(qdd) - C.dot(qd) - K.dot(q)


# lCols = ['HydroFxi_[N]','HydroFyi_[N]','HydroFzi_[N]','HydroMxi_[N-m]','HydroMyi_[N-m]','HydroMzi_[N-m]']
# dfPH = pd.DataFrame(data=np.column_stack((time,fh)), columns=['Time_[s]']+lCols)

# --- Plot Loads
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(time, fh_HD[:,5] , '-', label='HD')
ax.plot(time, fh_ED[:,5] , '--', label='ED')
ax.plot(time, fh_HD2[:,5], ':', label='HD2')
ax.set_xlabel('')
ax.set_ylabel('')
ax.legend()


# --- Plot
fig,axes = plt.subplots(6, 3, sharey=False, figsize=(12.8,8.5)) # (6.4,4.8)
fig.subplots_adjust(left=0.11, right=0.95, top=0.95, bottom=0.07, hspace=0.40, wspace=0.22)
# DOF
for iCol, (col1,col2) in enumerate(zip(qCol1, qCol2)):
    # Positions
    axes[iCol,0].plot(time, qPtfm_lin[:,iCol]  , '-',   label='Ptfm (lin)', lw=3, alpha=0.4)
    axes[iCol,0].plot(time, dfOF[col1].values  , '-', label='Hydro Ref')
    axes[iCol,0].plot(time, dfOF[col2].values  , '--', label='Ptfm')
    axes[iCol,0].plot(time, qPtfm[:,iCol]      , ':',   label='Ptfm (Transfered)')
    axes[iCol,0].plot(time, qRef [:,iCol]      , '-.',  label='Hydro ref (Transfered)')
    axes[iCol,0].set_ylabel(col1.replace('_',' '))
for iCol, (col1,col2) in enumerate(zip(qdCol1, qdCol2)):
    # Velocities
    axes[iCol,1].plot(time, qdPtfm_lin[:,iCol] , '-',   label='Ptfm (lin)', lw=3, alpha=0.4)
    axes[iCol,1].plot(time, dfOF[col1].values  , '-', label='Hydro Ref')
    axes[iCol,1].plot(time, dfOF[col2].values  , '--', label='Ptfm')
    axes[iCol,1].plot(time, qdPtfm[:,iCol]      , ':',   label='Ptfm (Transfered)')
    axes[iCol,1].plot(time, qdRef [:,iCol]      , '-.',  label='Hydro ref (Transfered)')
    axes[iCol,1].set_ylabel(col1.replace('_',' '))

for iCol, (col1,col2) in enumerate(zip(qddCol1, qddCol2)):
    # Acceleration
    axes[iCol,2].plot(time, qddPtfm_lin[:,iCol]  , '-',   label='Ptfm (lin)', lw=3, alpha=0.4)
    axes[iCol,2].plot(time, dfOF[col1].values  , '-', label='Hydro Ref')
    axes[iCol,2].plot(time, dfOF[col2].values  , '--', label='Ptfm')
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
