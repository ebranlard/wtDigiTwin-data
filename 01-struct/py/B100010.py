"""
Equations of motion
model name: B100010
"""
import numpy as np
from numpy import cos, sin, pi, sqrt
def info():
    """ Return information about current model present in this package """
    I=dict()
    I['name']='B100010'
    I['nq']=2
    I['nu']=0
    I['sq']=['x','phi_y']
    I['su']=[]
    return I

def forcing(t,q=None,qd=None,p=None,u=None,z=None):
    """ Non linear mass forcing 
    q:  degrees of freedom, array-like: ['x(t)', 'phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['M_B', 'g', 'x_BG', 'z_BG']
    u:  inputs, dictionary with keys: []
           where each values is a function of time
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    FF = np.zeros((2,1))
    FF[0,0] = p['M_B']*p['x_BG']*qd[1]**2*cos(q[1])+p['M_B']*p['z_BG']*qd[1]**2*sin(q[1])
    FF[1,0] = p['M_B']*p['g']*p['x_BG']*cos(q[1])+p['M_B']*p['g']*p['z_BG']*sin(q[1])
    return FF

def mass_matrix(q=None,p=None,z=None):
    """ Non linear mass matrix 
     q:  degrees of freedom, array-like: ['x(t)', 'phi_y(t)']
     p:  parameters, dictionary with keys: ['J_yy_B', 'M_B', 'x_BG', 'z_BG']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((2,2))
    MM[0,0] = p['M_B']
    MM[0,1] = -p['M_B']*p['x_BG']*sin(q[1])+p['M_B']*p['z_BG']*cos(q[1])
    MM[1,0] = -p['M_B']*p['x_BG']*sin(q[1])+p['M_B']*p['z_BG']*cos(q[1])
    MM[1,1] = p['J_yy_B']
    return MM

def M_lin(q=None,p=None,z=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'phi_y(t)']
    p:  parameters, dictionary with keys: ['J_yy_B', 'M_B', 'x_BG', 'z_BG']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((2,2))
    MM[0,0] = p['M_B']
    MM[0,1] = p['M_B']*(-p['x_BG']*sin(q[1])+p['z_BG']*cos(q[1]))
    MM[1,0] = p['M_B']*(-p['x_BG']*sin(q[1])+p['z_BG']*cos(q[1]))
    MM[1,1] = p['J_yy_B']
    return MM

def C_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear damping matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'phi_y(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: ['M_B', 'x_BG', 'z_BG']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    CC = np.zeros((2,2))
    CC[0,0] = 0
    CC[0,1] = -2*p['M_B']*p['x_BG']*qd[1]*cos(q[1])-2*p['M_B']*p['z_BG']*qd[1]*sin(q[1])
    CC[1,0] = 0
    CC[1,1] = 0
    return CC

def K_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear stiffness matrix 
    q:  degrees of freedom, array-like: ['x(t)', 'phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['M_B', 'g', 'x_BG', 'z_BG']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    KK = np.zeros((2,2))
    KK[0,0] = 0
    KK[0,1] = p['M_B']*p['x_BG']*qd[1]**2*sin(q[1])-p['M_B']*p['z_BG']*qd[1]**2*cos(q[1])
    KK[1,0] = 0
    KK[1,1] = p['M_B']*p['g']*p['x_BG']*sin(q[1])-p['M_B']*p['g']*p['z_BG']*cos(q[1])
    return KK

def B_lin(q=None,qd=None,p=None,u=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'phi_y(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: []
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    The columns of B correspond to:   []\\ 
    """
    BB = np.zeros((0,0))
    return BB

