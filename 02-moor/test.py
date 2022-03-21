import unittest
import os
import numpy as np    
import matplotlib.pyplot as plt
import welib
from welib.yams.models.simulator import *
from welib.yams.models.generator_oneRigidBody import generateOneRigidBodyModel
from welib.yams.models.generator import generateModel
from welib.tools.stats import mean_rel_err
from welib.FEM.utils import rigidTransformationTwoPoints, rigidTransformationTwoPoints_Loads
np.set_printoptions(linewidth=300, precision=5)

MyDir=os.path.dirname(__file__)
import platform


qop=None
qdop=None

packageDir  = os.path.join(MyDir,'py')

def YAMSsim(fstFilename, modelName, CG_on_z=False, qop=None, qdop=None, tMax=None, create=True, Refz=0, K_Moor0=None):
    fstFilename = os.path.join(MyDir,fstFilename)

    # Transform K_Moor from 0 point to refz
    P_HDRef = np.array((0,0,0))
    P_EDRef = np.array((0,0,Refz))
    T_HD2ED= rigidTransformationTwoPoints(P_HDRef, P_EDRef)
    T_ED2HD= rigidTransformationTwoPoints(P_EDRef, P_HDRef)
    T_HD2ED_l= rigidTransformationTwoPoints_Loads(P_HDRef, P_EDRef)
    T_ED2HD_l= rigidTransformationTwoPoints_Loads(P_EDRef, P_HDRef)
    K_Moor =T_HD2ED_l.dot(K_Moor0.dot(T_ED2HD))

    if create:
        if modelName[0]=='B':
            generateOneRigidBodyModel(modelName, CG_on_z=CG_on_z, packageDir=packageDir)
        else:
            generateModel(modelName, aero_forces=False, moor_loads=False, hydro_loads=False, packageDir=packageDir)

    WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)  # TODO
    sim=SimulatorFromOF(WT, modelName=modelName, packageDir=packageDir)
    if modelName[0]=='B':
        time, dfFS, p = sim.setupSim(tMax=tMax, flavor='onebody', J_at_Origin=True)
    else:
        time, dfFS, p = sim.setupSim(tMax=tMax, J_at_Origin=True)

    sim.qop=qop
    sim.qdop=qdop
    for i in range(6):
        for j in range(6):
            if j>=i:
                sim.p['KM_{}{}'.format(i,j)] = K_Moor[i,j]
    sim.p['z_T0'] = -Refz


    sim.simulate(out=False, prefix='')
    p=sim.p
    return sim

def comp(df1,df2,k):
    return mean_rel_err(y1=df1[k], y2=df2[k], method='meanabs', verbose=True, varname=k)

class Test(unittest.TestCase):

    def test_B000010_moorO(self):
        fstFilename='SparNoRNA_MD1HD0SD0_F000010/Main.fst'; modelName='B000010_moorO'; 
        Refz=25
        create = True

        K_Moor0=np.array(
        [[ 2.00000e+07, -0.00000e+00,  0.00000e+00, -0.00000e+00, -2.00000e+08,  0.00000e+00],
         [ 0.00000e+00,  0.00000e+00,  0.00000e+00,  2.00000e-01,  0.00000e+00,  0.00000e+00],
         [ 0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00],
         [ 0.00000e+00,  2.22000e+01,  0.00000e+00,  2.22200e+02,  0.00000e+00,  0.00000e+00],
         [-1.99997e+08,  0.00000e+00,  0.00000e+00,  0.00000e+00,  1.99987e+09, -0.00000e+00],
         [ 0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00,  0.00000e+00]])

        sim = YAMSsim(fstFilename, modelName, qop=qop, qdop=qdop, create=create, tMax=10, K_Moor0=K_Moor0, Refz=Refz)
        dfNL=sim.dfNL; dfLI=sim.dfLI; dfFS=sim.dfFS
        #sim.plot()
        # NL
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        # Lin
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);

    def test_B100010_moorO(self):
        fstFilename='Spar_MD1HD0SD0_F100010/Main.fst'; modelName='B100010_moorO'; 
        Refz=20
        create = True
        K_Moor0= np.array(
            [[ 6.53172e+05,-0.00000e+00,-1.25000e+01,-0.00000e+00, 5.80472e+06,-0.00000e+00],
             [ 0.00000e+00, 6.53154e+05, 0.00000e+00,-5.80442e+06, 0.00000e+00,-6.60000e+01],
             [-1.20000e+01,-0.00000e+00, 6.02505e+05, 0.00000e+00,-1.39400e+02,-0.00000e+00],
             [-0.00000e+00,-5.80970e+06,-0.00000e+00, 5.15665e+08, 9.97000e+01,-9.96950e+03],
             [ 5.81000e+06,-0.00000e+00,-2.34460e+03, 0.00000e+00, 5.15673e+08,-0.00000e+00],
             [-0.00000e+00,-6.55000e+01,-0.00000e+00, 2.32050e+03, 2.32000e+01, 4.96495e+08]])
        sim = YAMSsim(fstFilename, modelName, qop=qop, qdop=qdop, create=create, tMax=30, K_Moor0=K_Moor0, Refz=Refz)
        dfNL=sim.dfNL; dfLI=sim.dfLI; dfFS=sim.dfFS

        sim.plot()
        # NL
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 19.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 19.0);
        # Lin
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 3.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 3.0);


    def test_F5T1N0S1(self):
        fstFilename='Spar_MD1HD0SD0_F5T1N0S1/Main.fst'; modelName='F5T1N0S1_fnd'; 
        qop=[0,0.00,0,-0.0*np.pi/180, 0.0*np.pi/180,0,0]; aero_forces=False;  qdop=[0,0,0,0,0,0,10/60*2*np.pi];
        Refz=20

        create=True
        if platform.node()=='ebranlar-36947s':
            create=False

        K_Moor0=np.array(
        [[ 1.77778e+04, -0.00000e+00, -1.00000e-01, -0.00000e+00,  9.27501e+04, -0.00000e+00],
         [-0.00000e+00,  1.77776e+04, -0.00000e+00, -9.27462e+04, -0.00000e+00, -1.60000e+00],
         [-1.00000e-01,  0.00000e+00,  1.51606e+04, -0.00000e+00, -1.20000e+00,  0.00000e+00],
         [-0.00000e+00, -9.27772e+04,  0.00000e+00,  3.64207e+07,  1.20000e+00, -1.22200e+02],
         [ 9.27811e+04,  0.00000e+00, -1.53000e+01, -0.00000e+00,  3.64210e+07,  0.00000e+00],
         [-0.00000e+00, -1.60000e+00, -0.00000e+00,  5.16000e+01,  5.00000e-01,  3.44370e+07]])
        sim = YAMSsim(fstFilename, modelName, qop=qop, qdop=qdop, create=create, tMax=10, K_Moor0=K_Moor0, Refz=Refz)

        dfNL=sim.dfNL; dfLI=sim.dfLI; dfFS=sim.dfFS
        sim.plot()
        # NL
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 2.0);
        k='PtfmSway_[m]';    eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 2.0);
        k='PtfmRoll_[deg]';  eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 2.0);
        k='PtfmYaw_[deg]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 20.0);
        k='Azimuth_[deg]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='RotSpeed_[rpm]';  eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='Q_TFA1_[m]';      eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 7.0);
        # Lin
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 3.0);
        k='PtfmSway_[m]';    eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 2.0);
        k='PtfmRoll_[deg]';  eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 3.0);
        k='PtfmYaw_[deg]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 11.0);
        k='Azimuth_[deg]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='RotSpeed_[rpm]';  eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='Q_TFA1_[m]';      eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 12.0);

if __name__ == '__main__':
    Test().test_B000010_moorO()
    Test().test_B100010_moorO()
    #Test().test_F5T1N0S1()
    #unittest.main()
    plt.show()
