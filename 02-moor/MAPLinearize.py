from welib.moor.mappp import Map
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(linewidth=300, precision=4)



gravity=9.80665
WtrDens=1025
WtrDepth=220

# gravity=0.01
# WtrDens=0.1


moor = Map('../TetraSparModel/MAP.dat', WtrDepth=WtrDepth, gravity=gravity, WtrDens=WtrDens)

# --- Plot initial 
fig, ax = moor.plot(numPoints=20, colors=['k'], ls='--')

# --- Linearization with no displacement (NOTE: is it in equilibrium?)
epsilon = 1e-3 # finite difference epsilon
K = moor.linear(epsilon)    
print("\nLinearized stiffness matrix with 0.0 vessel displacement:\n")
print(np.around(K,2))
 
# --- Linearization with a given surge displacement
surge = 0.0 # 5 meter surge displacements
heave = 0.0 # 5 meter surge displacements
moor.displace_vessel(surge,0,heave,0,0,0)
moor.update_states(t=0.0, dt=0)
K = moor.linear(epsilon)    
print("\nLinearized stiffness matrix with %2.2f surge vessel displacement:\n"%(surge))
print(K)

plt.show()
