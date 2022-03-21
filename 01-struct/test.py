import unittest
import os
import numpy as np    
import matplotlib.pyplot as plt

from welib.yams.models.simulator import *
from welib.yams.models.generator_oneRigidBody import generateOneRigidBodyModel
from welib.yams.models.generator import generateModel
from welib.tools.stats import mean_rel_err
np.set_printoptions(linewidth=300, precision=5)

MyDir=os.path.dirname(__file__)
import platform
create=True
if platform.node()=='ebranlar-36947s':
    create=True
    create=False

def YAMSsim(fstFilename, modelName, CG_on_z=False, qop=None, qdop=None, tMax=None, create=True):
    fstFilename = os.path.join(MyDir,fstFilename)
    packageDir  = os.path.join(MyDir,'py')

    if create:
        if modelName[0]=='B':
            generateOneRigidBodyModel(modelName, CG_on_z=CG_on_z)
        else:
            generateModel(modelName, aero_forces=False, moor_loads=False, hydro_loads=False)

    WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)  # TODO


    sim=SimulatorFromOF(WT, modelName=modelName, packageDir=packageDir)
    if modelName[0]=='B':
        time, dfFS, p = sim.setupSim(tMax=tMax, flavor='onebody', J_at_Origin=True)
    else:
        time, dfFS, p = sim.setupSim(tMax=tMax, J_at_Origin=True)

    sim.qop=qop
    sim.qdop=qdop
    sim.simulate(out=False, prefix='')
    p=sim.p
    return sim

def comp(df1,df2,k):
    return mean_rel_err(y1=df1[k], y2=df2[k], method='meanabs', verbose=True, varname=k)

class Test(unittest.TestCase):

    def test_F5T1N0S1(self):

        fstFilename = 'Spar_F5T1N0S1/Main_Spar_ED.fst'; modelName='F5T1N0S1_fnd'; 
        qop=[0,0,0,0,0,0,0]; qdop=[0,0,0,0,0,0,10/60*2*np.pi];
        sim = YAMSsim(fstFilename, modelName, qop=qop, qdop=qdop, create=create, tMax=10)

        dfNL=sim.dfNL
        dfLI=sim.dfLI
        dfFS=sim.dfFS

        # NL
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 5.0);
        k='PtfmSway_[m]';    eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 3.0);
        k='PtfmRoll_[deg]';  eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 7.0);
        k='PtfmYaw_[deg]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 2.0);
        k='Azimuth_[deg]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='RotSpeed_[rpm]';  eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        # Lin
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 3.0);
        k='PtfmSway_[m]';    eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 5.0);
        k='PtfmRoll_[deg]';  eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 5.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 4.0);
        k='PtfmYaw_[deg]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='Azimuth_[deg]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='RotSpeed_[rpm]';  eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);

        sim.plot(export=False)
        plt.show()

    def test_F5T1RNA(self):

        fstFilename = 'Spar_F5T1RNA/Main_Spar_ED.fst'; modelName='F5T1RNA_fnd'; 
        qop=[0,0,0,0,0,0]; qdop=[0,0,0,0,0,0];
        sim = YAMSsim(fstFilename, modelName, qop=qop, qdop=qdop, create=create, tMax=10)

        dfNL=sim.dfNL
        dfLI=sim.dfLI
        dfFS=sim.dfFS

        # NL
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 5.0);
        k='PtfmSway_[m]';    eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmRoll_[deg]';  eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 4.0);
        k='PtfmYaw_[deg]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 2.0);
        # Lin
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 5.0);
        k='PtfmSway_[m]';    eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='PtfmRoll_[deg]';  eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 2.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 5.0);
        k='PtfmYaw_[deg]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);

        sim.plot(export=False)
        plt.show()


if __name__ == '__main__':
    Test().test_F5T1RNA()
    #unittest.main()
