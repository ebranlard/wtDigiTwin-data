"""
Equations of motion
model name: B101010_hydroO
"""
import numpy as np
from numpy import cos, sin, pi, sqrt
def info():
    """ Return information about current model present in this package """
    I=dict()
    I['name']='B101010_hydroO'
    I['nq']=3
    I['nu']=3
    I['sq']=['x','z','phi_y']
    I['su']=['F_hx','F_hz','M_hy']
    return I

def forcing(t,q=None,qd=None,p=None,u=None,z=None):
    """ Non linear mass forcing 
    q:  degrees of freedom, array-like: ['x(t)', 'z(t)', 'phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['M_B', 'g', 'x_BG', 'z_BG']
    u:  inputs, dictionary with keys: ['F_hx', 'F_hz', 'M_hy']
           where each values is a function of time
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    FF = np.zeros((3,1))
    FF[0,0] = p['M_B']*p['x_BG']*qd[2]**2*cos(q[2])+p['M_B']*p['z_BG']*qd[2]**2*sin(q[2])+u['F_hx'](t,q,qd)
    FF[1,0] = -p['M_B']*p['g']-p['M_B']*p['x_BG']*qd[2]**2*sin(q[2])+p['M_B']*p['z_BG']*qd[2]**2*cos(q[2])+u['F_hz'](t,q,qd)
    FF[2,0] = p['M_B']*p['g']*p['x_BG']*cos(q[2])+p['M_B']*p['g']*p['z_BG']*sin(q[2])+u['M_hy'](t,q,qd)
    return FF

def mass_matrix(q=None,p=None,z=None):
    """ Non linear mass matrix 
     q:  degrees of freedom, array-like: ['x(t)', 'z(t)', 'phi_y(t)']
     p:  parameters, dictionary with keys: ['J_yy_B', 'M_B', 'x_BG', 'z_BG']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((3,3))
    MM[0,0] = p['M_B']
    MM[0,1] = 0
    MM[0,2] = -p['M_B']*p['x_BG']*sin(q[2])+p['M_B']*p['z_BG']*cos(q[2])
    MM[1,0] = 0
    MM[1,1] = p['M_B']
    MM[1,2] = -p['M_B']*p['x_BG']*cos(q[2])-p['M_B']*p['z_BG']*sin(q[2])
    MM[2,0] = -p['M_B']*p['x_BG']*sin(q[2])+p['M_B']*p['z_BG']*cos(q[2])
    MM[2,1] = -p['M_B']*p['x_BG']*cos(q[2])-p['M_B']*p['z_BG']*sin(q[2])
    MM[2,2] = p['J_yy_B']
    return MM

def M_lin(q=None,p=None,z=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'z(t)', 'phi_y(t)']
    p:  parameters, dictionary with keys: ['J_yy_B', 'M_B', 'x_BG', 'z_BG']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((3,3))
    MM[0,0] = p['M_B']
    MM[0,1] = 0
    MM[0,2] = p['M_B']*(-p['x_BG']*sin(q[2])+p['z_BG']*cos(q[2]))
    MM[1,0] = 0
    MM[1,1] = p['M_B']
    MM[1,2] = p['M_B']*(-p['x_BG']*cos(q[2])-p['z_BG']*sin(q[2]))
    MM[2,0] = p['M_B']*(-p['x_BG']*sin(q[2])+p['z_BG']*cos(q[2]))
    MM[2,1] = p['M_B']*(-p['x_BG']*cos(q[2])-p['z_BG']*sin(q[2]))
    MM[2,2] = p['J_yy_B']
    return MM

def C_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear damping matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'z(t)', 'phi_y(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: ['M_B', 'x_BG', 'z_BG']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    CC = np.zeros((3,3))
    CC[0,0] = 0
    CC[0,1] = 0
    CC[0,2] = -2*p['M_B']*p['x_BG']*qd[2]*cos(q[2])-2*p['M_B']*p['z_BG']*qd[2]*sin(q[2])
    CC[1,0] = 0
    CC[1,1] = 0
    CC[1,2] = 2*p['M_B']*p['x_BG']*qd[2]*sin(q[2])-2*p['M_B']*p['z_BG']*qd[2]*cos(q[2])
    CC[2,0] = 0
    CC[2,1] = 0
    CC[2,2] = 0
    return CC

def K_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear stiffness matrix 
    q:  degrees of freedom, array-like: ['x(t)', 'z(t)', 'phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['M_B', 'g', 'x_BG', 'z_BG']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    KK = np.zeros((3,3))
    KK[0,0] = 0
    KK[0,1] = 0
    KK[0,2] = p['M_B']*p['x_BG']*qd[2]**2*sin(q[2])-p['M_B']*p['z_BG']*qd[2]**2*cos(q[2])
    KK[1,0] = 0
    KK[1,1] = 0
    KK[1,2] = p['M_B']*p['x_BG']*qd[2]**2*cos(q[2])+p['M_B']*p['z_BG']*qd[2]**2*sin(q[2])
    KK[2,0] = 0
    KK[2,1] = 0
    KK[2,2] = p['M_B']*p['g']*p['x_BG']*sin(q[2])-p['M_B']*p['g']*p['z_BG']*cos(q[2])
    return KK

def B_lin(q=None,qd=None,p=None,u=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'z(t)', 'phi_y(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: []
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    The columns of B correspond to:   [F_hx(t), F_hz(t), M_hy(t)]\\ 
    """
    BB = np.zeros((3,3))
    BB[0,0] = 1
    BB[0,1] = 0
    BB[0,2] = 0
    BB[1,0] = 0
    BB[1,1] = 1
    BB[1,2] = 0
    BB[2,0] = 0
    BB[2,1] = 0
    BB[2,2] = 1
    return BB

