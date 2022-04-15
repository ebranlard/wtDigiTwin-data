import unittest
import os
import numpy as np    
import matplotlib.pyplot as plt
import welib
from welib.yams.models.simulator import *
from welib.yams.models.generator_oneRigidBody import generateOneRigidBodyModel
from welib.yams.models.generator import generateModel
from welib.tools.stats import mean_rel_err
from welib.fast.hydrodyn_driver import hydroSimFromOpenFAST, hydroSimLinFromOpenFAST
from welib.fast.hydrodyn import HydroDyn
import pickle
np.set_printoptions(linewidth=300, precision=5)

MyDir=os.path.dirname(__file__)
packageDir  = os.path.join(MyDir,'py')

import platform


def comp(df1,df2,k):
    return mean_rel_err(y1=df1[k], y2=df2[k], method='meanabs', verbose=True, varname=k)

def prescribedLoadSim(fstFilename, modelName, create, tMax=None):
    fstFilename = os.path.join(MyDir, fstFilename)
    packageDir = os.path.join(MyDir, 'py')
    if create:
        if modelName[0]=='B':
            generateOneRigidBodyModel(modelName, packageDir=packageDir)
        else:
            generateModel(modelName, aero_forces=False, moor_loads=False, hydro_loads=True, packageDir=packageDir)

    # --- Run non linear and linear simulation using a FAST model as input
    print('----------------------- SETUP SIMULATION -----------------------------------------')
    WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)  # TODO
    sim = SimulatorFromOF(WT, modelName=modelName, packageDir=packageDir)
    if modelName[0]=='B':
        time, dfFS, p = sim.setupSim(tMax=tMax, flavor='onebody', J_at_Origin=True)
    else:
        time, dfFS, p = sim.setupSim(tMax=tMax, J_at_Origin=True)

    sim.setPrescribedHydroInputs()
    #uop=None
    #du=None
    qop  = None
    qdop = None
    # Using mean as op
    # qop  = np.array([np.mean(dfFS[c]) for c in WT.q_channels])
    # qdop = np.array([np.mean(dfFS[c]) for c in WT.qd_channels])*0
    # --- Simulation
    #sim.setInputs(u, du, uop, qop, qdop)
    sim.simulate(out=False)
    return sim


class TestHydro(unittest.TestCase):

    def test_PrescribeMotion_SparF101010(self):

        fstFilename = os.path.join(MyDir, 'SparNoRNA_MD0HD1SD0_F101010/Main.fst')

        tMax=10
        #tMax=None

        q0= np.zeros(6)

        # --- Generate linearized model
        hd = HydroDyn(fstFilename)
        MCKF0 = hd.linearize_RigidMotion2Loads(q0=q0)
        # --- Time integration (TODO, use simulator?)
        print('Time integration linear hydrodynamics ')
        dfLI, _, _    =    hydroSimLinFromOpenFAST(fstFilename, tMax=tMax, MCKF=MCKF0, q0=q0)
        print('Time integration nonlinear hydrodynamics ')
        dfNL, dfOF, msy, _ = hydroSimFromOpenFAST(fstFilename, tMax=tMax, verbose=False)
        # NL
        k='HydroFxi_[N]';   eps=comp(dfNL, dfOF, k); self.assertLessEqual(eps, 6.0);
        k='HydroFzi_[N]';   eps=comp(dfNL, dfOF, k); self.assertLessEqual(eps, 1.0);
        k='HydroMyi_[N-m]'; eps=comp(dfNL, dfOF, k); self.assertLessEqual(eps, 1.0);
        print('Lin')

        k='HydroFxi_[N]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 6.0);
        k='HydroFzi_[N]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='HydroMyi_[N-m]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 7.0);

    def test_PrescribeMotion_SparF111111(self):

        #fstFilename = 'SparNoRNA_F101010/Main.fst'; # Takes 6min - Should give really good match (0 signas have noise)
        fstFilename = os.path.join(MyDir, 'Spar_MD0HD1SD0_F111111/Main.fst');

        tMax=10
        #tMax=None

        q0= np.zeros(6)

        # --- Generate linearized model
        hd = HydroDyn(fstFilename)
        MCKF0 = hd.linearize_RigidMotion2Loads(q0=q0)
        # --- Time integration (TODO, use simulator?)
        print('Time integration linear hydrodynamics ')
        dfLI, _ ,_    =    hydroSimLinFromOpenFAST(fstFilename, tMax=tMax, MCKF=MCKF0, q0=q0)
        print('Time integration nonlinear hydrodynamics ')
        dfNL, dfOF, msy, _ = hydroSimFromOpenFAST(fstFilename, tMax=tMax, verbose=False)
        # NL
        k='HydroFxi_[N]';   eps=comp(dfNL, dfOF, k); self.assertLessEqual(eps, 1.0);
        k='HydroFyi_[N]';   eps=comp(dfNL, dfOF, k); self.assertLessEqual(eps, 1.0);
        k='HydroFzi_[N]';   eps=comp(dfNL, dfOF, k); self.assertLessEqual(eps, 1.0);
        k='HydroMxi_[N-m]'; eps=comp(dfNL, dfOF, k); self.assertLessEqual(eps, 1.0);
        k='HydroMyi_[N-m]'; eps=comp(dfNL, dfOF, k); self.assertLessEqual(eps, 1.0);
        k='HydroMzi_[N-m]'; eps=comp(dfNL, dfOF, k); self.assertLessEqual(eps, 1.0);
        print('Lin')

        k='HydroFzi_[N]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='HydroMxi_[N-m]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 10.0);
        k='HydroMyi_[N-m]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 10.0);


    # --------------------------------------------------------------------------------}
    # --- Prescribed loads 
    # --------------------------------------------------------------------------------{
    def test_PrescribeLoads_SparF101010(self):
        tMax=10
        #modelName = 'B111111_hydro0'; fstFilename = 'Spar_MD0HD1SD0_F111111/Main.fst'; hydro=True
        # modelName = 'B111111_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F111111/Main.fst'; hydro=True
        modelName = 'B101010_hydro0'; fstFilename = 'SparNoRNA_MD0HD1SD0_F101010/Main.fst'; 
        sim = prescribedLoadSim(fstFilename, modelName, create=False, tMax=tMax)
        dfNL=sim.dfNL; dfLI=sim.dfLI; dfFS=sim.dfFS
        # sim.plot()
        # plt.show()
        # NL
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmHeave_[m]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        # Lin
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='PtfmHeave_[m]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 1.0);

    def test_PrescribeLoads_SparF111111(self):
        tMax=10
        modelName = 'B111111_hydro0'; fstFilename = 'Spar_MD0HD1SD0_F111111/Main.fst'; hydro=True
        sim = prescribedLoadSim(fstFilename, modelName, create=False, tMax=tMax)
        dfNL=sim.dfNL; dfLI=sim.dfLI; dfFS=sim.dfFS
        #sim.plot()
        # NL
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmSway_[m]';    eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmHeave_[m]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmRoll_[deg]';  eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 1.0);
        k='PtfmYaw_[deg]';   eps=comp(dfNL, dfFS, k); self.assertLessEqual(eps, 12.0);
        # Lin
        k='PtfmSurge_[m]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 5.0);
        k='PtfmSway_[m]';    eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 2.0);
        k='PtfmHeave_[m]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 5.0);
        k='PtfmRoll_[deg]';  eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 3.0);
        k='PtfmPitch_[deg]'; eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 4.0);
        k='PtfmYaw_[deg]';   eps=comp(dfNL, dfLI, k); self.assertLessEqual(eps, 9.0);
        #sim.plot()
        #plt.show()




if __name__ == '__main__':
    #TestHydro().test_PrescribeMotion_SparF101010()
    #TestHydro().test_PrescribeMotion_SparF111111()
    TestHydro().test_PrescribeLoads_SparF101010()
    TestHydro().test_PrescribeLoads_SparF111111()
    #unittest.main()
    plt.show()
