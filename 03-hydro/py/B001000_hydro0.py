"""
Equations of motion
model name: B001000_hydro0
"""
import numpy as np
from numpy import cos, sin, pi, sqrt
def info():
    """ Return information about current model present in this package """
    I=dict()
    I['name']='B001000_hydro0'
    I['nq']=1
    I['nu']=1
    I['sq']=['z']
    I['su']=['F_hz']
    return I

def forcing(t,q=None,qd=None,p=None,u=None,z=None):
    """ Non linear mass forcing 
    q:  degrees of freedom, array-like: ['z(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['M_B', 'g']
    u:  inputs, dictionary with keys: ['F_hz']
           where each values is a function of time
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    FF = np.zeros((1,1))
    FF[0,0] = -p['M_B']*p['g']+u['F_hz'](t,q,qd)
    return FF

def mass_matrix(q=None,p=None,z=None):
    """ Non linear mass matrix 
     q:  degrees of freedom, array-like: ['z(t)']
     p:  parameters, dictionary with keys: ['M_B']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((1,1))
    MM[0,0] = p['M_B']
    return MM

def M_lin(q=None,p=None,z=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['z(t)']
    p:  parameters, dictionary with keys: ['M_B']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((1,1))
    MM[0,0] = p['M_B']
    return MM

def C_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear damping matrix 
    q:  degrees of freedom at operating point, array-like: ['z(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: []
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    CC = np.zeros((1,1))
    CC[0,0] = 0
    return CC

def K_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear stiffness matrix 
    q:  degrees of freedom, array-like: ['z(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: []
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    KK = np.zeros((1,1))
    KK[0,0] = 0
    return KK

def B_lin(q=None,qd=None,p=None,u=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['z(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: []
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    The columns of B correspond to:   [F_hz(t)]\\ 
    """
    BB = np.zeros((1,1))
    BB[0,0] = 1
    return BB

