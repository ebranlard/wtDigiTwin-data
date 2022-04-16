"""  
 Pendulum modelling a pinned rigid body under gravity and hydro loads

 Everything is expressed at the "reference" point "0" compared to body origin.
 NOTE: HydroDyn return forces at (0,0,0), so "z_B0" is the distance from the body origin to the global (0,0,0) point
       But if bodies moves, these forces are still at the origin...
 
 Equation of motions:
  
     JO theta_ddot  = M *g *z_OG sin\theta  + Mhy  + z_B0 (F_hx\cos\theta - F_hz\sin\theta)    with z_OG<0
    


"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.system.mech_system import MechSystem
import welib.weio as weio
from welib.tools.clean_exceptions import *
MyDir=os.path.dirname(__file__)
from welib.yams.windturbine import FASTWindTurbine
from welib.FEM.utils import transferRigidLoads
from welib.yams.utils import transferLoadsZPoint

# --- Script parameters

useYAMSparams=True

noHydro=True
noHydro=False

loadsAtO = False



if noHydro:
    fstFilename = '../01-struct/SparNoRNA_F000010_NoHydro/Main.fst'
else:
    fstFilename = 'SparNoRNA_MD0HD1SD0_F000010_NoCdCpCa/Main.fst'
    fstFilename =  'TS_MD0HD1SD0_F000010_NoCdCpCa/Main.fst';



# --- Open Main channels from OpenFAST simulations
dfFS = weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
time =dfFS['Time_[s]'].values
if noHydro:
    vMy   = time*0
    vF_hx = time *0
    vF_hz = time *0
else:
    vMy   = dfFS['HydroMyi_[N-m]']
    vF_hx = dfFS['HydroFxi_[N]']
    vF_hz = dfFS['HydroFzi_[N]']
vtheta     =  dfFS['Q_P_[rad]'].values # [rad] solution
vphi_x     =  dfFS['Q_R_[rad]'].values # [rad] solution
vphi_y     =  dfFS['Q_P_[rad]'].values # [rad] solution


# --- Parameters
if useYAMSparams:
    WT = FASTWindTurbine(fstFilename, algo = 'OpenFAST')
    p  = WT.yams_parameters(flavor='onebody', J_at_Origin=True)
    Mtot = p['M_B']
    g    = p['g']
    JO   = p['J_yy_B']
    x_OG = p['x_BG']
    z_OG = p['z_BG']
    z_O0 = p['z_B0']
else:
    Mtot = 6.530370E+06 # [kg]     body mass
    g    = 10  # [m/s-2]  acceleration of gravity
    z_OG = -25 -1.038637E+02
    z_O0 = -25
    JO = 1.193506E+11
    x_OG =0

theta0    =  vtheta[0]     # Initial angle    [rad]
thetadot0 =  0*np.pi/180 # Initial velocity [rad/s]

# Matrices
M=np.array([[JO]])
K=np.array([[0]])
C=np.array([[0]])



if loadsAtO:
    cols = ['HydroFxi_[N]', 'HydroFyi_[N]', 'HydroFzi_[N]', 'HydroMxi_[N-m]', 'HydroMyi_[N-m]', 'HydroMzi_[N-m]']
    zRef =  p['z_B0'] # TODO TODO TODO WEIRD SHOULD BE -
    print('zRef',zRef)
    P_HDRef = np.array((0,0,0))
    P_EDRef = np.array((0,0,zRef))
    lM = dfFS[cols].values
    #MT = transferRigidLoads(lM.T, P_HDRef, P_EDRef).T
    MT = transferLoadsZPoint(lM.T, zRef, vphi_x, vphi_y).T

    dfFS2 = pd.DataFrame(data=MT, columns=cols)
    vMy2  = dfFS2['HydroMyi_[N-m]'].values
    vMy3  = vMy + z_O0 * (vF_hx*np.cos(vtheta) - vF_hz*np.sin(vtheta))
    #vMy3  = vMy + z_O0 * vF_hx

#     fig,axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
#     fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
#     ax = axes[0]
#     ax.plot(time, vMy  , '-'  , label='vMy')
#     ax.plot(time, vMy2 , '--' , label='vMy transfer')
#     ax.plot(time, vMy3 , ':'  , label='vMy non linear')
#     ax = axes[1]
# #     ax.plot(time, vF_hz  , '-'  , label='vFz')
# #     ax.plot(time, MT[:,2] , '--' , label='vFz transfer')
# #     ax.plot(time, vF_hx   , '-'  , label='vFx')
# #     ax.plot(time, MT[:,0] , '--' , label='vFx transfer')
# #     ax.plot(time, dfFS ['HydroMxi_[N-m]'] , '-'  , label='vMx')
# #     ax.plot(time, dfFS2['HydroMxi_[N-m]'] , '--' , label='vMx transfer')
#     ax.plot(time, dfFS ['HydroMyi_[N-m]'] , '-'  , label='vMy')
#     ax.plot(time, dfFS2['HydroMyi_[N-m]'] , '--' , label='vMy transfer')
# #     ax.plot(time, dfFS ['HydroMzi_[N-m]'] , '-'  , label='vMy')
# #     ax.plot(time, dfFS2['HydroMzi_[N-m]'] , '--' , label='vMy transfer')
#     ax.set_xlabel('')
#     ax.set_ylabel('')
#     ax.legend()
#     plt.show()





# --- Definition of RHS
def Fx(t,x,xdot):
    """ Returns force on pendulum"""
    theta      = x[0,0]
    Force      = np.zeros((1,1))
    if loadsAtO:
        Myh2 = np.interp(t, time, vMy2)
        Myh = Myh2
    else:
        Myh  = np.interp(t, time, vMy)
        F_hx = np.interp(t, time, vF_hx)
        F_hz = np.interp(t, time, vF_hz)
        Myh1 = z_O0 * (F_hx*np.cos(theta) - F_hz*np.sin(theta))
        #Myh1 = z_O0 * F_hx
        Myh += Myh1
    Force[0,0] = g* Mtot*(x_OG*np.cos(theta) + z_OG*np.sin(theta))  + Myh
    #import pdb; pdb.set_trace()
    return Force

# --- Define a system and perform time integration
sys=MechSystem(M, C, K, Fx, x0=[theta0], xdot0=[thetadot0] )
res=sys.integrate(time, method='RK45') # **options):

fig, axes = sys.plot()

axes[0].plot(time, dfFS['Q_P_[rad]']  , 'k:')
axes[1].plot(time, dfFS['QD_P_[rad/s]'], 'k:')

if __name__=="__main__":
    plt.show()
if __name__=="__test__":
    pass

