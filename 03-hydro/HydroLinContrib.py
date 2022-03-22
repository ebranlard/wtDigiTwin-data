""" 
Linearize python Hydrodyn module wrt motion at ref point
# NOTE: ballast and MG have zero contrib
# HydroI: has zero contrib
# Buoyancy: contribute to K
# HydroA: contribute to M
# HydroD: contribute to C
# End:    Mostly K, but contributes to Ch, Mh in term 3,3 (z,z) (from rest position) but more 3:5, 3:5 if tilted
# 
# M: Mostly HydroA
# C: Mostly HydroD
# K: Buoyancy and End
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pickle
# Local 
from welib.weio.fast_input_file import FASTInputFile
from welib.fast.hydrodyn import HydroDyn
from welib.fast.fast_mesh import *
from welib.tools.clean_exceptions import *
from welib.tools.tictoc import *
from welib.system.linearization import numerical_jacobian
from welib.tools.colors import fColrs

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import style
# style.use('ggplot')


def main(fstFilename, base='TetraSpar', q0=np.zeros(6)):
    # --- Initialize a python HydroDyn instance
    hd = HydroDyn(fstFilename)

    # --- Linearization
    qd0  = np.zeros(6)
    qdd0 = np.zeros(6)
    dq   = [0.01]*3 + [0.01]*3
    dqd  = [0.01]*3 + [0.01]*3
    dqdd = [0.1]*3 + [0.1]*3

    pdefault ={'Buoyancy':False, 'HydroA':False, 'HydroD':False, 'HydroI':False, 'End':False}

    splitContrib = {}

    Msum = np.zeros((6,6))
    Csum = np.zeros((6,6))
    Ksum = np.zeros((6,6))
    for k in pdefault.keys():
        p=pdefault.copy()
        p[k]=True
        print('---------------------------------------  ',k)
        Mh, Ch, Kh = hd.linearize_RigidMotion2Loads(q0, optsM=p)
        print(np.around(Kh,4))
        print(np.around(Ch,4))
        print(np.around(Mh,4))
        splitContrib[k] = (Mh, Ch, Kh)
        Msum+= Mh
        Csum+= Ch
        Ksum+= Kh

    with Timer('Linearize all contrib'):
        Mall, Call, Kall = hd.linearize_RigidMotion2Loads(q0)
    allContrib = (Mall, Call, Kall)
    print('--------------------------------------- ALL CONTRIB')
    print(np.around(Kall,4))
    print(np.around(Call,4))
    print(np.around(Mall,4))
    print('--------------------------------------- SUM CONTRIB')
    print(np.around(Ksum,4))
    print(np.around(Csum,4))
    print(np.around(Msum,4))
    print('--------------------------------------- DIFF')
    print(Kall-Ksum)
    print(Call-Csum)
    print(Mall-Msum)
#     pickle.dump((allContrib, splitContrib),  open('Hlin_{}.pkl'.format(base),'wb'))
# 
    allContrib, splitContrib = pickle.load(open('Hlin_{}.pkl'.format(base),'rb'))

    comps=['M','C','K']
    iC=2
    base0=base
    for scaled in [True,False]:
        if scaled:
            base=base0+'_scaled'
        else:
            base=base0
        for iC in [0,1,2]:
            fig,axes = plt.subplots(6, 6, sharey=False, figsize=(11.4,11.4)) # (6.4,4.8)
            fig.subplots_adjust(left=0.16, right=0.95, top=0.93, bottom=0.05, hspace=0.5, wspace=0.5)
            max1=np.max(np.abs(allContrib[iC][:3,:]))*1.05
            max2=np.max(np.abs(allContrib[iC][3:,:]))*1.05
            for i in np.arange(6):
                for j in np.arange(6):
                    axes[i,j].bar(0, allContrib[iC][i,j], label='all', color='k')
                    bottom=0
                    for ii,(k,v) in enumerate(splitContrib.items()):
        #                 axes[i,j].bar(1,    v[iC][i,j], label=k, bottom=bottom, color=fColrs(ii))
                        axes[i,j].bar(ii+2, v[iC][i,j], color=fColrs(ii), label=k)
                        bottom+=v[iC][i,j]
                    if scaled:
                        if i<3:
                            axes[i,j].set_ylim([-max1, max1])
                        else:
                            axes[i,j].set_ylim([-max2, max2])
                    axes[i,j].set_xlabel('')
                    axes[i,j].set_ylabel('')
                    axes[i,j].set_title('{},{}'.format(i+1,j+1), fontsize=11)
                    if i==0 and j==0:
                        axes[i,j].legend(bbox_to_anchor=(-0.12, 0.9))
                fig.suptitle(comps[iC]+base, fontsize=16)
                fig.savefig('_figs/{}_{}.png'.format(base,comps[iC]))
    # 




np.set_printoptions(linewidth=300, precision=4)

q0   = np.zeros(6)
main('../TetraSparModel/TetraSpar_SWT-3p6-130_RigidPtfm.fst', base='TetraSpar', q0=q0)
# 
q0   = [1,1,1,0.1,0.1,0.1]
main('../TetraSparModel/TetraSpar_SWT-3p6-130_RigidPtfm.fst', base='TetraSpar_titled', q0=q0)


plt.show()



if __name__ == '__main__':
    pass
