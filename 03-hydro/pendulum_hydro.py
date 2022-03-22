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


noHydro=True
noHydro=False

if noHydro:
    fstFilename = os.path.join(MyDir, '../models/SparNoRNA_F000010_NoHydro/Main.fst')
else:
    fstFilename = os.path.join(MyDir, '../models/SparNoRNA_F000010/Main.fst')

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

# --- Parameters
Mtot = 6.530370E+06 # [kg]     body mass
g    = 10  # [m/s-2]  acceleration of gravity
z_OG = -25 -1.038637E+02
z_O0 = -25
theta0    =  dfFS['PtfmPitch_[deg]'].values[0]*np.pi/180 # Initial angle    [rad]
thetadot0 =  0*np.pi/180 # Initial velocity [rad/s]
JO = 1.193506E+11

# Matrices
M=np.array([[JO]])
K=np.array([[0]])
C=np.array([[0]])

# --- Definition of RHS
def Fx(t,x,xdot):
    """ Returns force on pendulum"""
    Myh = np.interp(t, time, vMy)
    F_hx = np.interp(t, time, vF_hx)
    F_hz = np.interp(t, time, vF_hz)
    theta      = x[0,0]
    Force      = np.zeros((1,1))
    Myh1 = z_O0 * (F_hx*np.cos(theta) - F_hz*np.sin(theta))
    #Myh1 = z_O0 * F_hx
    Myh += Myh1
    Force[0,0] = g* Mtot*z_OG*np.sin(theta)  + Myh
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

