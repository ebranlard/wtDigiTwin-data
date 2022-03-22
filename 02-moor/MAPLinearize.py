from welib.moor.mappp import Map
import matplotlib.pyplot as plt
import numpy as np
from welib.FEM.utils import rigidTransformationTwoPoints, rigidTransformationTwoPoints_Loads
np.set_printoptions(linewidth=300, precision=4)



gravity=9.80665
WtrDens=1025
WtrDepth=220

Refz=16
# gravity=0.01
# WtrDens=0.1


moor = Map('../TetraSparModel/MAP.dat', WtrDepth=WtrDepth, gravity=gravity, WtrDens=WtrDens)

# --- Plot initial 
fig, ax = moor.plot(numPoints=20, colors=['k'], ls='-')

# --- Linearization with no displacement (NOTE: is it in equilibrium?)
epsilon = 1e-3 # finite difference epsilon
K = moor.linear(epsilon)    
print("\nLinearized stiffness matrix at (0,0,0):\n")
print(np.around(K,1))

P_HDRef = np.array((0,0,0))
P_EDRef = np.array((0,0,Refz))
T_ED2HD  = rigidTransformationTwoPoints(P_EDRef, P_HDRef)
T_HD2ED_l= rigidTransformationTwoPoints_Loads(P_HDRef, P_EDRef)
K_P = T_HD2ED_l.dot(K.dot(T_ED2HD))

print("\nLinearized stiffness matrix at (0,0,RefZ):\n")
print(np.around(K_P,1))


 
# --- Linearization with a given surge displacement
surge = 0.0 # 5 meter surge displacements
heave = 0.0 # 5 meter surge displacements
moor.displace_vessel(surge,0,heave,0,0,0)
moor.update_states(t=0.0, dt=0)
K = moor.linear(epsilon)    
print("\nLinearized stiffness matrix with %2.2f surge vessel displacement:\n"%(surge))
print(K)

plt.show()
