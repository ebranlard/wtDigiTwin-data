"""  
 Heave pendulum

 Everything is expressed at the "reference" point 
 
 Equation of motions:
  
     (M+Ma) z_ddot  = -M *g + Fz
    


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

model='SparNoRNA'
model='TS'
# model='Spar'

if model=='SparNoRNA':
    fstFilename = os.path.join(MyDir, 'SparNoRNA_MD0HD1SD0_F001000/Main.fst')
elif model=='Spar':
    fstFilename = os.path.join(MyDir, 'Spar_MD0HD1SD0_F001000/Main.fst')
elif model=='TS':
    fstFilename = os.path.join(MyDir, 'TS_MD0HD1SD0_F001000/Main.fst')
    fstFilename = os.path.join(MyDir, 'TS_MD0HD1SD0_F001000_NoCdCpCa/Main.fst')
dfFS  = weio.read(fstFilename.replace('.fst','.outb')).toDataFrame()
time  = dfFS['Time_[s]'].values
vF_hz = dfFS['HydroFzi_[N]']

# --- Parameters
if model=='SparNoRNA':
    # Mtot = 6.530370E+06 # [kg]     body mass
    # Ma   = 0.1e6
    Ma =0
    M = 6.518904e+06
    M = 6.530e+06
    g    = 10  # [m/s-2]  acceleration of gravity
    z0    = -1 # Initial heave    [m]
    zdot0 =  0 # Initial velocity [m/s]
elif model=='Spar':
    Ma =0
    M = 8067146.9397829715
    g    = 10  # [m/s-2]  acceleration of gravity
    z0    = -0.1 # Initial heave    [m]
    zdot0 =  0 # Initial velocity [m/s]
elif model=='TS':
    g    = 9.80665 # [m/s-2]  acceleration of gravity
    M = 5584174.89141556
    #M  = 5580970
    Ma = 0
    z0    = -0.1 # Initial heave    [m]
    zdot0 =  0   # Initial velocity [m/s]

# z_OG = -25 -1.038637E+02
# z_O0 = -25

# Matrices
M=np.array([[M+Ma]])
K=np.array([[0]])
C=np.array([[0]])

# --- Definition of RHS
def Fx(t,x,xdot):
    """ Returns force on pendulum"""
    F_hz = np.interp(t, time, vF_hz)
    z        = x[0,0]
    Force      = np.zeros((1,1))
    Force[0,0] = -g* M  + F_hz
    return Force

# --- Define a system and perform time integration
sys=MechSystem(M, C, K, Fx, x0=[z0], xdot0=[zdot0] )
res=sys.integrate(time, method='RK45') # **options):

fig, axes = sys.plot()

axes[0].plot(time, dfFS['Q_Hv_[m]']  , 'k:', label='OpenFAST')
axes[1].plot(time, dfFS['QD_Hv_[m/s]'], 'k:')
axes[0].legend()

if __name__=="__main__":
    plt.show()
if __name__=="__test__":
    pass

