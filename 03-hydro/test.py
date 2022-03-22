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


class TestHydro(unittest.TestCase):

    def test_PrescribeMotion_SparF101010(self):

        fstFilename = os.path.join(MyDir, 'SparNoRNA_F101010/Main.fst')

        tMax=10
        #tMax=None

        q0= np.zeros(6)

        # --- Generate linearized model
        hd = HydroDyn(fstFilename)
        MCKF0 = hd.linearize_RigidMotion2Loads(q0=q0)
        # --- Time integration (TODO, use simulator?)
        print('Time integration linear hydrodynamics ')
        dfLI, _    =    hydroSimLinFromOpenFAST(fstFilename, tMax=tMax, MCK=MCKF0[:3], q0=q0, F0=MCKF0[-1])
        print('Time integration nonlinear hydrodynamics ')
        dfNL, dfOF, msy = hydroSimFromOpenFAST(fstFilename, tMax=tMax, verbose=False)
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
        dfLI, _    =    hydroSimLinFromOpenFAST(fstFilename, tMax=tMax, MCK=MCKF0[:3], q0=q0, F0=MCKF0[-1])
        print('Time integration nonlinear hydrodynamics ')
        dfNL, dfOF, msy = hydroSimFromOpenFAST(fstFilename, tMax=tMax, verbose=False)
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


if __name__ == '__main__':
    #TestHydro().test_PrescribeMotion_SparF101010()
    #TestHydro().test_PrescribeMotion_SparF111111()
    unittest.main()
    plt.show()
