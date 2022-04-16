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

    # --- Spar No RNA
    #fstFilename = 'SparNoRNA_F101010/Main.fst'; # Takes 6min - Should give really good match (0 signas have noise)
    fstFilename = 'SparNoRNA_MD0HD1SD0_F000010/Main.fst';
#     fstFilename = 'SparNoRNA_MD0HD1SD0_F000010_NoRef/Main.fst';
#     fstFilename = 'SparNoRNA_MD0HD1SD0_F000010_NoRef_NoMor/Main.fst';
#     fstFilename = 'SparNoRNA_MD0HD1SD0_F111111/Main.fst'; # NL:OK, Lin:pretty good

    # --- Spar
#     fstFilename = 'Spar_MD0HD1SD0_F111111/Main.fst'; # NL:OK, 
    fstFilename = 'Spar_MD0HD1SD0_F111111_NoCd/Main.fst'; # NL:OK, Lin:Pseudo-OK (e.g. yaw)

    # --- Tetra Spar
#     fstFilename = 'TS_MD0HD1SD0_F001000/Main.fst'; 
#     fstFilename = 'TS_MD0HD1SD0_F101010/Main.fst'; # NL:Perfcet, Lin: No ( - Should give perfect match (NOTE: No Vn*2/3)
#     fstFilename = 'TS_MD0HD1SD0_F101010_NoCdCpCa/Main.fst'; 
#     fstFilename = 'TS_MD0HD1SD0_F101010_NoCdCp/Main.fst';
#     fstFilename = 'TS_MD0HD1SD0_F101010_NoCd/Main.fst';
#     fstFilename = 'TS_MD0HD1SD0_F111111/Main.fst'; # NL:Perfect, Lin: wrong becasue of Cd (NOTE: No Vn*2/3)
# #     fstFilename = 'TS_MD0HD1SD0_F111111_NoCdCpCa/Main.fst';
# #     fstFilename = 'TS_MD0HD1SD0_F111111_NoCdCp/Main.fst';
#     fstFilename = 'TS_MD0HD1SD0_F111111_NoCd/Main.fst'; # NL:OK, Lin:pretty good (e.g. init of Mzi is -5000 instead of -4000)

    tMax=1
    tMax=5
#     tMax=20
    tMax=None

    q0= np.zeros(6)

    # --- Generate linearized model
    print('>>>> FST',fstFilename)
    hd = HydroDyn(fstFilename)
    zRef=0 # <<<<<<<<<<<<<<<<<<<<<<< TODO TODO TODO TODO TODO HACK
    MCKF0 = hd.linearize_RigidMotion2Loads(q0=q0, RefPointMotion=(0,0,zRef), RefPointMapping=(0,0,zRef))
    pickle.dump(MCKF0, open(fstFilename.replace('.fst','_hydroPyPrescrMotion_Lin.pkl'),'wb'))
    # OR read it
    MCKF0 = pickle.load(open(fstFilename.replace('.fst','_hydroPyPrescrMotion_Lin.pkl'),'rb'))

    print('Mh\n',np.around(MCKF0[0]        , 0))
    print('Ch\n',np.around(MCKF0[1]        , 0))
    print('Kh\n',np.around(MCKF0[2]        , 0))
    print('F0'  ,np.around(MCKF0[3].ravel(), 0))

    fig=None

    print('Time integration nonlinear hydrodynamics ')
    dfPH, dfOF, msy, fig = hydroSimFromOpenFAST(fstFilename, tMax=tMax, verbose=False, fig=fig)
    print('Time integration linear hydrodynamics ')
    dfLI, _, fig    =    hydroSimLinFromOpenFAST(fstFilename, tMax=tMax, MCKF=MCKF0, q0=q0, fig=fig)


    plt.show()

