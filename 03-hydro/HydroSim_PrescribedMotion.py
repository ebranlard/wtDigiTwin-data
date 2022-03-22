""" 
Compute hydrodynamics loads on a structure under rigid body motion using "Python HydroDyn"

- Hydrodynamic parameters are read from an HydroDyn input file
- Motion (displacement, velocities, accelerations) are taken from an OpenFAST simulation
  NOTES:
      - requires latest dev branch of OpenFAST to get platform reference outputs for now
      - if ascii files are used, the output resolution "OutFmt" should have sufficient digits)
- The Motion is applied at the HydroDyn reference point 
- Loads are computed at the HydroDyn reference point using the python hydrodyn module.

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# Local 
import welib.weio as weio
from welib.fast.hydrodyn_driver import hydroSimFromOpenFAST, hydroSimLinFromOpenFAST
from welib.fast.hydrodyn import HydroDyn
import pickle

MyDir=os.path.dirname(__file__)
np.set_printoptions(linewidth=300, precision=2)

if __name__ == '__main__':
    #fstFilename = 'SparNoRNA_F101010/Main.fst'; # Takes 6min - Should give really good match (0 signas have noise)
    #fstFilename = 'TS_MD0HD1SD0_F101010/Main.fst'; # Takes 7min - Should give perfect match (NOTE: No Vn*2/3)
    fstFilename = 'TS_MD0HD1SD0_F111111/Main.fst'; # Takes 7min - Should give perfect match (NOTE: No Vn*2/3)

    tMax=10
    tMax=None

    q0= np.zeros(6)

    # --- Generate linearized model
    hd = HydroDyn(fstFilename)
    MCKF0 = hd.linearize_RigidMotion2Loads(q0=q0)
    pickle.dump(MCKF0, open(fstFilename.replace('.fst','_hydroPyPrescrMotion_Lin.pkl'),'wb'))
    # OR read it
    MCKF0 = pickle.load(open(fstFilename.replace('.fst','_hydroPyPrescrMotion_Lin.pkl'),'rb'))

    print('Mh\n',np.around(MCKF0[0]        , 0))
    print('Ch\n',np.around(MCKF0[1]        , 0))
    print('Kh\n',np.around(MCKF0[2]        , 0))
    print('F0'  ,np.around(MCKF0[3].ravel(), 0))

    print('Time integration linear hydrodynamics ')
    dfLI, _    =    hydroSimLinFromOpenFAST(fstFilename, tMax=tMax, MCK=MCKF0[:3], q0=q0, F0=MCKF0[-1])
    print('Time integration nonlinear hydrodynamics ')
    dfPH, dfOF, msy = hydroSimFromOpenFAST(fstFilename, tMax=tMax, verbose=False)


    plt.show()

# if __name__ == '__test__':
#     fstFilename = os.path.join(MyDir, '../../../data/Spar/Main_Spar_ED_HydroExample.fst');
#     dfPH, dfOF = hydroSim(fstFilename, tMax=0.1, plot=False, json=False)
#     np.testing.assert_almost_equal(dfPH['HydroFxi_[N]'].values  /1e6,dfOF['HydroFxi_[N]'].values  /1e6, 5)
#     np.testing.assert_almost_equal(dfPH['HydroFyi_[N]'].values  /1e6,dfOF['HydroFyi_[N]'].values  /1e6, 5)
#     np.testing.assert_almost_equal(dfPH['HydroFzi_[N]'].values  /1e6,dfOF['HydroFzi_[N]'].values  /1e6, 5)
#     np.testing.assert_almost_equal(dfPH['HydroMxi_[N-m]'].values/1e6,dfOF['HydroMxi_[N-m]'].values/1e6, 3)
#     np.testing.assert_almost_equal(dfPH['HydroMyi_[N-m]'].values/1e6,dfOF['HydroMyi_[N-m]'].values/1e6, 3)
#     np.testing.assert_almost_equal(dfPH['HydroMzi_[N-m]'].values/1e6,dfOF['HydroMzi_[N-m]'].values/1e6, 3)
