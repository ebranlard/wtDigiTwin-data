import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
# Local 

from numpy.linalg import inv
from welib.FEM.utils import rigidTransformationTwoPoints18, rigidTransformationOnePointToPoints18
from welib.FEM.utils import rigidTransformationTwoPoints, rigidTransformationTwoPoints_Loads
from welib.yams.rotations import BodyXYZ_A, SmallRot_DCM
from welib.fast.postpro import find_matching_pattern
from welib.fast.hydrodyn import HydroDyn
from welib.fast.hydrodyn_driver import *
from welib.fast.linmodel import matToSIunits, subMat, matLabelNoUnit, matLabelReplace
from welib.tools import *

import weio

def cleanMat(df):
    df = matLabelNoUnit (df)
    df = matLabelReplace(df,'Ptfm','')
    df = matLabelReplace(df,'HDMorison','')
    return df

np.set_printoptions(linewidth=300, precision=2, threshold=10000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


fstFilename='./HydroLin/TetraSparRigid_Lin_Hydro.fst'; i=1; Read=False; Refz=16
fstFilename='./SparNoRNA_F11111_MiniStatesOutputs/Main.fst'; i=1; Read=True; Refz=25
fstFilename='./SparNoRNA_F11111_MiniStatesOutputs/MainNoRef.fst'; i=1; Read=True; Refz=0


# --- Python model of HydroDyn
hd = HydroDyn(fstFilename)

# --- Read lin file
if Read:
    with Timer('Read'):
        lin = weio.read(fstFilename.replace('.fst','.{}.lin'.format(i)))
    pickle.dump(lin, open(fstFilename.replace('.fst','.1.pkl'), 'wb'))
lin = pickle.load(open(fstFilename.replace('.fst','.1.pkl'), 'rb'))

with Timer('toDataFrame'):
    dfs = lin.toDataFrame()

#return self['x_info']['Description'])
# print('----------- X ')
# [print(sx) for sx in lin.xdescr()]
# print('----------- Y ')
# [print(sy) for sy in lin.ydescr()]
# print('----------- U ')
# [print(su) for su in lin.udescr()]
# print(len(lin.udescr()))
# print(len(lin['u']))
# print(lin.keys())
# print(lin['x'])
# print(dfs.keys())
# print(lin['dUdy'].shape)

_, nodes = find_matching_pattern(lin.udescr(), 'HDMorisonTxN(\d+)_\[m\]')
nNodes = np.max(nodes)


# Inputs: HydroDyn/Morison nodes motion
uMorNds = []
for n in range(1,nNodes+1):
    uMorNds+=['HDMorisonTxN{}_[m]'.format(n), 'HDMorisonTyN{}_[m]'.format(n), 'HDMorisonTzN{}_[m]'.format(n), 'HDMorisonRxN{}_[rad]'.format(n), 'HDMorisonRyN{}_[rad]'.format(n), 'HDMorisonRzN{}_[rad]'.format(n), 'HDMorisonTVxN{}_[m/s]'.format(n), 'HDMorisonTVyN{}_[m/s]'.format(n), 'HDMorisonTVzN{}_[m/s]'.format(n), 'HDMorisonRVxN{}_[rad/s]'.format(n), 'HDMorisonRVyN{}_[rad/s]'.format(n), 'HDMorisonRVzN{}_[rad/s]'.format(n), 'HDMorisonTAxN{}_[m/s^2]'.format(n), 'HDMorisonTAyN{}_[m/s^2]'.format(n), 'HDMorisonTAzN{}_[m/s^2]'.format(n), 'HDMorisonRAxN{}_[rad/s^2]'.format(n), 'HDMorisonRAyN{}_[rad/s^2]'.format(n), 'HDMorisonRAzN{}_[rad/s^2]'.format(n)]

# Inputs: HydroDyn Ref Point motion
uHDRefP=['HDPtfm-RefPtTxN1_[m]', 'HDPtfm-RefPtTyN1_[m]', 'HDPtfm-RefPtTzN1_[m]', 'HDPtfm-RefPtRxN1_[rad]', 'HDPtfm-RefPtRyN1_[rad]', 'HDPtfm-RefPtRzN1_[rad]', 'HDPtfm-RefPtTVxN1_[m/s]', 'HDPtfm-RefPtTVyN1_[m/s]', 'HDPtfm-RefPtTVzN1_[m/s]', 'HDPtfm-RefPtRVxN1_[rad/s]', 'HDPtfm-RefPtRVyN1_[rad/s]', 'HDPtfm-RefPtRVzN1_[rad/s]', 'HDPtfm-RefPtTAxN1_[m/s^2]', 'HDPtfm-RefPtTAyN1_[m/s^2]', 'HDPtfm-RefPtTAzN1_[m/s^2]', 'HDPtfm-RefPtRAxN1_[rad/s^2]', 'HDPtfm-RefPtRAyN1_[rad/s^2]', 'HDPtfm-RefPtRAzN1_[rad/s^2]']

# Outputs: ElastoDyn Ptfm motion
yEDRefP=['PtfmTxN1_[m]', 'PtfmTyN1_[m]', 'PtfmTzN1_[m]', 'PtfmRxN1_[rad]', 'PtfmRyN1_[rad]', 'PtfmRzN1_[rad]', 'PtfmTVxN1_[m/s]', 'PtfmTVyN1_[m/s]', 'PtfmTVzN1_[m/s]', 'PtfmRVxN1_[rad/s]', 'PtfmRVyN1_[rad/s]', 'PtfmRVzN1_[rad/s]', 'PtfmTAxN1_[m/s^2]', 'PtfmTAyN1_[m/s^2]', 'PtfmTAzN1_[m/s^2]', 'PtfmRAxN1_[rad/s^2]', 'PtfmRAyN1_[rad/s^2]', 'PtfmRAzN1_[rad/s^2]']

# Outputs: hydro loads
yHydroF=['HDHydroFxi_[N]', 'HDHydroFyi_[N]', 'HDHydroFzi_[N]', 'HDHydroMxi_[Nm]', 'HDHydroMyi_[Nm]', 'HDHydroMzi_[Nm]']


#--- Extract the relevant part of the Jacobian
# print(list(dfs['dUdy'].index[-4:])) # Inputs
# print(list(dfs['dUdy'].columns[-4:]))
dUdy =  -cleanMat(subMat(dfs['dUdy'], uMorNds, yEDRefP)) # 18n x 18

# --- Rigid body transformation matrix from ED PtfmRef to HD nodes
P_EDRef = np.array((0,0,Refz))
P_HDRef = np.array((0,0,0))
T_ED2HDNodes = rigidTransformationOnePointToPoints18(P_EDRef, hd.u['Morison']['Mesh'].Position) # Positions are wrt to HDRef
print(hd.u['Morison']['Mesh'].Position)
T_HD2ED = rigidTransformationTwoPoints(P_HDRef, P_EDRef)

i=137;
print('T{}\n'.format(i+1),dUdy[i*18:(i+1)*18].round())
print('T{}\n'.format(i+1),T_ED2HDNodes[i*18:(i+1)*18].round())
# print('dUdy\n',dUdy.round())
# TT=dUdy.copy()
# # TT.values=T
# print('T   \n',TT.round())
# print('T   \n',T-dUdy.values)
print('T   \n',np.max(np.abs(T_ED2HDNodes-dUdy.values)))
# print(dUdy)
# print(dUdy.shape)


#--- Extract the relevant part of the Jacobian
dUdu = subMat(dfs['dUdu'], uMorNds, uHDRefP) # 18n x 18


# import pdb; pdb.set_trace()
# --- Extract the relevant part of the D matrix
Dtilde = cleanMat(subMat(dfs['D'], yHydroF, uMorNds))  # 6 x 18n

D=Dtilde.values.dot(dUdy.values)
# D=Dtilde.dot(dUdu)
# print(D.round())
D=Dtilde.dot(dUdy)
print('--- D matrix (K C M')
print(D.round())
Kh = D.iloc[:,0:6]
Ch = D.iloc[:,6:12]
Mh = D.iloc[:,12:18]
print('---- K, C, M ED to Hydro')
print(np.round(Kh))
print(np.round(Ch))
print(np.round(Mh))
print('---- K, C, M HD to Hydro')
print(np.round(Kh.dot(T_HD2ED)))
print(np.round(Ch.dot(T_HD2ED)))
print(np.round(Mh.dot(T_HD2ED)))

print('----------------')
# T_HD2ED = rigidTransformationTwoPoints_Loads(P_HDRef, P_EDRef)
print(np.round(Kh))
print(np.round(Kh.dot(T_HD2ED)))
print(T_HD2ED)


# raise


# --------------------------------------------------------------------------------}
# --- PYTHON 
# --------------------------------------------------------------------------------{
# --- Initialize a python HydroDyn instance
dfOF = weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
if 'PRPSurge_[m]' in dfOF.columns:
    qCol   = ['PRPSurge_[m]'    ,'PRPSway_[m]'    ,'PRPHeave_[m]'   ,'PRPRoll_[rad]'    ,'PRPPitch_[rad]'   ,'PRPYaw_[rad]'     ]
q0=dfOF[qCol].values[0,:]
lCols = ['HydroFxi_[N]','HydroFyi_[N]','HydroFzi_[N]','HydroMxi_[N-m]','HydroMyi_[N-m]','HydroMzi_[N-m]']
print('q0',q0)



# --- Linearization
with Timer('Linearize all contrib'):
    Mall, Call, Kall, f0 = hd.linearize_RigidMotion2Loads(q0=q0)

# --- py Lin at HydroRef Point
# base ='TetraSpar'
# allContrib, splitContrib = pickle.load(open('../../code/Hlin_{}.pkl'.format(base),'rb'))
# Mall, Call, Kall =allContrib
print('f0',f0)
print('f0',dfOF[lCols].values[0,:])

print(np.round(-Kall))
print(np.round(-Call))
print(np.round(-Mall))
# 
# with Timer('Non linear sim'):
#     dfLI, _         = hydroSimLinFromOpenFAST(fstFilename, tMax=None, MCK=(Mall,Call,Kall), q0=q0, F0=f0)
# with Timer('Non linear sim'):
#     dfPH, dfOF, msy = hydroSimFromOpenFAST(fstFilename, tMax=None, verbose=False)
# 
# 
# plt.show()


# # --- Convert the units to SI units
# #print('------- Convert to SI units')
# Minv=matToSIunits(Minv, 'D')
# #print(Minv)
# 
# # ---  Inverse the matrix
# #print('------- Inverse the matrix')
# M=inv(Minv)
# print(M)
# 
# 
# 
# # Rigid transformation matrix between DOFs of node j and k where node j is the leader node (Pk-Pj)
# Refz=16
# P1 = np.array((0,0,0))
# P2 = np.array((0,0,Refz))
# T= rigidTransformationTwoPoints(P2, P1)
# print(T)
# 
# 
# if __name__ == '__main__':
#     pass
