"""
Equations of motion
model name: B000010_hydro0
"""
import numpy as np
from numpy import cos, sin, pi, sqrt
def info():
    """ Return information about current model present in this package """
    I=dict()
    I['name']='B000010_hydro0'
    I['nq']=1
    I['nu']=3
    I['sq']=['phi_y']
    I['su']=['F_hx','F_hz','M_hy']
    return I

def forcing(t,q=None,qd=None,p=None,u=None,z=None):
    """ Non linear mass forcing 
    q:  degrees of freedom, array-like: ['phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['M_B', 'g', 'z_B0', 'z_BG']
    u:  inputs, dictionary with keys: ['F_hx', 'F_hz', 'M_hy']
           where each values is a function of time
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    FF = np.zeros((1,1))
    FF[0,0] = p['M_B']*p['g']*p['z_BG']*sin(q[0])+p['z_B0']*u['F_hx'](t,q,qd)*cos(q[0])-p['z_B0']*u['F_hz'](t,q,qd)*sin(q[0])+u['M_hy'](t,q,qd)
    return FF

def mass_matrix(q=None,p=None,z=None):
    """ Non linear mass matrix 
     q:  degrees of freedom, array-like: ['phi_y(t)']
     p:  parameters, dictionary with keys: ['J_yy_B']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((1,1))
    MM[0,0] = p['J_yy_B']
    return MM

def M_lin(q=None,p=None,z=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['phi_y(t)']
    p:  parameters, dictionary with keys: ['J_yy_B']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((1,1))
    MM[0,0] = p['J_yy_B']
    return MM

def C_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear damping matrix 
    q:  degrees of freedom at operating point, array-like: ['phi_y(t)']
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
    q:  degrees of freedom, array-like: ['phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['M_B', 'g', 'z_B0', 'z_BG']
    u:  inputs at operating point, dictionary with keys: ['F_hx', 'F_hz']
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    KK = np.zeros((1,1))
    KK[0,0] = -p['M_B']*p['g']*p['z_BG']*cos(q[0])-p['z_B0']*(-u['F_hx']*sin(q[0])-u['F_hz']*cos(q[0]))
    return KK

def B_lin(q=None,qd=None,p=None,u=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['phi_y(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: ['z_B0']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    The columns of B correspond to:   [F_hx(t), F_hz(t), M_hy(t)]\\ 
    """
    BB = np.zeros((1,3))
    BB[0,0] = p['z_B0']*cos(q[0])
    BB[0,1] = -p['z_B0']*sin(q[0])
    BB[0,2] = 1
    return BB

