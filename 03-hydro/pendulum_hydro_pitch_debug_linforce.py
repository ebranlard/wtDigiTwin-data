import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from numpy import cos, sin
# Local 
import weio
import py.B000010_hydro0 as model
from welib.yams.windturbine import FASTWindTurbine
from welib.tools.clean_exceptions import *
from welib.system.mech_system import MechSystem
MyDir=os.path.dirname(__file__)

# --- OpenFAST
fstFilename = os.path.join(MyDir, 'SparNoRNA_MD0HD1SD0_F000010/Main.fst')
if os.path.exists(fstFilename.replace('.fst','.outb')): 
    dfFS = weio.read(fstFilename.replace('.fst','.outb')).toDataFrame()
    time =dfFS['Time_[s]'].values
elif os.path.exists(fstFilename.replace('.fst','.out')): 
    dfFS = weio.read(fstFilename.replace('.fst','.out')).toDataFrame()
    time =dfFS['Time_[s]'].values

F_hx        = dfFS['HydroFxi_[N]'   ].values
F_hz        = dfFS['HydroFzi_[N]'   ].values
M_hy        = dfFS['HydroMyi_[N-m]' ].values
theta_y     = dfFS['PtfmPitch_[deg]'].values*np.pi/180
theta_y_dot = dfFS['QD_P_[rad/s]'].values

#
WT = FASTWindTurbine(fstFilename, twrShapes=[0,2], nSpanTwr=50)
p = WT.yams_parameters(flavor='onebody', J_at_Origin=True)

print(WT)

# --- Non linear model
# Initial conditions
q0  = np.array([theta_y[0]])
qd0 = np.array([theta_y_dot[0]])
# Inputs
u=dict()
u['F_hx'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFxi_[N]']  )
u['F_hz'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroFzi_[N]']  )
u['M_hy'] = lambda t,q=None,qd=None: np.interp(t, time, dfFS['HydroMyi_[N-m]'])

# --- Linear model
qop  = np.array([0])
qop  = np.array([np.mean(theta_y)])
qdop = np.array([0])
uop  = {'F_hx':np.mean(F_hx), 'F_hz':np.mean(F_hz), 'M_hy':np.mean(M_hy)}
dq0  = q0  - qop
dqd0 = qd0 - qdop

M_lin   = model.M_lin(qop,p)
C_lin   = model.C_lin(qop,qdop,p,uop)
K_lin   = model.K_lin(qop,qdop,p,uop) 
B_lin   = model.B_lin(qop,qdop,p,uop)


t=0
MM0      = model.mass_matrix(q0,p)
forcing0 = model.forcing(t,q0,qd0,p,u)

# --- Integrate non-linear system
# fM = lambda x: model.mass_matrix(x, p)
# fF = lambda t,x,xd: model.forcing(t, x, xd, p=p, u=u)
# sysNL = MechSystem(fM, F=fF, x0=q0, xdot0=qd0 )
# resNL = sysNL.integrate(time, method='RK45')


FNL = np.zeros(len(time))
Fg = p['M_B']*p['g']*p['z_BG']*sin(theta_y)
Fx =                                        p['z_B0']*F_hx*cos(theta_y)
Fz =                                                                   -p['z_B0']*F_hz*sin(theta_y)
FM = M_hy
F = Fg+Fx+Fz+FM

# for t in time:
#     F = p['M_B']*p['g']*p['z_BG']*sin(q[0])+p['z_B0']*u['F_hx'](t,q,qd)*cos(q[0])-p['z_B0']*u['F_hz'](t,q,qd)*sin(q[0])+u['M_hy'](t,q,qd)


# --- Integrate linear system
# fdu = lambda t,x=None,xd=None: interpArray(t, time, du)
# fF  = lambda t,x=None,xd=None: B_lin.dot( interpArray(t, time, du) )
# 
# sysLI = MechSystem(M=M_lin, K=K_lin, C=C_lin, F=fF, x0=dq0, xdot0=dqd0)
# resLI=sysLI.integrate(time, method='RK45') # **options):
# sysLI._B = B_lin

print('B_lin',B_lin)

FKlin = -K_lin[0,0] * theta_y
Fxlin = B_lin[0,0] * F_hx
Fzlin = B_lin[0,1] * F_hz
FMlin = B_lin[0,2] * M_hy

Flin = FMlin + FKlin + Fxlin + Fzlin




fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(time, F        , label='Non linear')
ax.plot(time, Flin, '--', label='linear')
ax.plot(time, Flin+Fz, 'k:',label='linear+Fz')
# ax.plot(time, Fg   , label='Fg')
# ax.plot(time, FKlin, label='FK', ls='--')
# 
# ax.plot(time, Fx   , label='Fx    ')
# ax.plot(time, Fxlin, label='Fx lin', ls='--')
# 
# ax.plot(time, Fz   , label='Fz    ')
# ax.plot(time, Fzlin, label='Fz lin', ls='--')
# 
# ax.plot(time, FM   , label='My    ')
# ax.plot(time, FMlin, label='My lin', ls='--')

ax.set_xlabel('')
ax.set_ylabel('')
ax.legend()
plt.show()  
