# 
# Pinned uniform rigid body floater with hydrodynamics 
# 
# NOTE: due to hydro, the two cases (where the structure is shifted) won't give the same
# 
# NOTE: added hub mass and overhang!!!
#
#      The reason is starting at a given pitch angle, SubDyn ended up crashing around pitch=0!!! TO DEBUG!
#
# 
# Equation of motions:
#  
#     JO theta_ddot + M *g *l/2 sin\theta =   tau_Hydro
#
import numpy as np

# --- Parameters
Mtot = 5.4e6 # [kg]     floater mass
L    = 150   # [m]      floater length
rho  = 7850  # [kg/m^3] material density 
g    = 10    # [m/s-2] acceleration of gravity

# --- Derevied parameters
D  = np.sqrt(4*Mtot/(rho*np.pi*L))
A  = np.pi*D**2
JO = Mtot * L**2/3                 # moment of inertia wrt to O (around x or y axis) [kg m^2]
JG = Mtot * L**2/12                # moment of inertia wrt to G (around x or y axis) [kg m^2]
Jz = Mtot * (D/2)**2/2             # moment of inertia           around z-axis [kg m^2]

zG = -L/2

omega = np.sqrt(Mtot * g * L/2 / JO)


print('theta= ' ,10*np.pi/180)
print('D    = ' ,D)
print('t    = ' ,D/2)
print('A    = ' ,A)
print('m(z) = ' ,rho*A)
print('JOy=JOx = ' ,JO)
print('JGy=JGx = ' ,JG)
print('Jz      = ' ,Jz)
print('m z_G = ',Mtot*zG)
print('omega = ',omega )
print('f     = ',omega/(2*np.pi) )
print('T     = ',2*np.pi/omega )
