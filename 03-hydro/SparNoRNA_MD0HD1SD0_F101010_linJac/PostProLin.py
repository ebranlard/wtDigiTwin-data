"""


Map++:
  will return the stiff


"""







import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio

from welib.fast.linmodel import matToSIunits
from numpy.linalg import inv
from welib.FEM.utils import rigidTransformationTwoPoints
from welib.yams.rotations import BodyXYZ_A, SmallRot_DCM





fstFilename='./Main.fst'; i=1;


# --- Read lin file
lin = weio.read(fstFilename.replace('.fst','.{}.lin'.format(i)))
dfs = lin.toDataFrame()
#return self['x_info']['Description'])
print('----------- X ')
[print(sx) for sx in lin.xdescr()]
print('----------- Y ')
[print(sy) for sy in lin.ydescr()]
print('----------- U ')
[print(su) for su in lin.ydescr()]
print(lin.keys())
print(lin['x'])




dfD = dfs['D']
# --- Extract the relevant 6x6 matrix
#print('------- Extract 6x6 matrix')
# Cols  = ['HDHydroFxi_[N]', 'HDHydroFyi_[N]', 'HDHydroFzi_[N]', 'HDHydroMxi_[N-m]', 'HDHydroMyi_[N-m]', 'HDHydroMzi_[N-m]']
# Cols  = ['PtfmFxN1_[N]', 'PtfmFyN1_[N]', 'PtfmFzN1_[N]', 'PtfmMxN1_[Nm]', 'PtfmMyN1_[Nm]', 'PtfmMzN1_[Nm]']
# Lines = ['PtfmTAxt_[m/s^2]', 'PtfmTAyt_[m/s^2]', 'PtfmTAzt_[m/s^2]', 'PtfmRAxt_[deg/s^2]', 'PtfmRAyt_[deg/s^2]', 'PtfmRAzt_[deg/s^2]']
# Lines = ['PtfmTAxt_[m/s^2]', 'PtfmTAyt_[m/s^2]', 'PtfmTAzt_[m/s^2]', 'PtfmRAxt_[deg/s^2]', 'PtfmRAyt_[deg/s^2]', 'PtfmRAzt_[deg/s^2]']
# Cols  = ['PtfmFzN1_[N]']
# Lines = ['PtfmTAzt_[m/s^2]']
# 
# 
# missingRows = [l for l in Lines if l not in dfD.index]
# missingCols = [c for c in Cols  if c not in dfD.columns]
# if len(missingRows)>0:
#     raise Exception('The following rows are missing from outputs: {}'.format(missingRows))
# if len(missingCols)>0:
#     raise Exception('The following columns are missing from inputs: {}'.format(missingCols))
# Minv = dfD[Cols].loc[Lines].copy()
# #print('Minv\n',Minv)
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
# Rigid transformation matrix between DOFs of node j and k where node j is the leader node (Pk-Pj)
Refz=16
P1 = np.array((0,0,0))
P2 = np.array((0,0,Refz))
T= rigidTransformationTwoPoints(P2, P1)
print(T)
# 

if __name__ == '__main__':
    pass
