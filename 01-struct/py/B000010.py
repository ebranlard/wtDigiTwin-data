"""
Equations of motion
model name: B000010
"""
import numpy as np
from numpy import cos, sin, pi, sqrt
def info():
    """ Return information about current model present in this package """
    I=dict()
    I['name']='B000010'
    I['nq']=1
    I['nu']=4
    I['sq']=['phi_y']
    I['su']=['Derivative(omega_x, t)','Derivative(omega_z, t)','omega_x','omega_z']
    return I

def forcing(t,q=None,qd=None,p=None,u=None,z=None):
    """ Non linear mass forcing 
    q:  degrees of freedom, array-like: ['phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['J_xx_B', 'J_zz_B', 'M_B', 'g', 'x_BG', 'z_BG']
    u:  inputs, dictionary with keys: ['omega_x', 'omega_z']
           where each values is a function of time
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    FF = np.zeros((1,1))
    FF[0,0] = -p['J_xx_B']*u['omega_x'](t,q,qd)*u['omega_z'](t,q,qd)+p['J_zz_B']*u['omega_x'](t,q,qd)*u['omega_z'](t,q,qd)+p['M_B']*p['g']*p['x_BG']*cos(q[0])+p['M_B']*p['g']*p['z_BG']*sin(q[0])
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
    p:  parameters, dictionary with keys: ['M_B', 'x_BG', 'y_BG', 'z_BG']
    u:  inputs at operating point, dictionary with keys: ['omega_x', 'omega_z']
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    CC = np.zeros((1,1))
    CC[0,0] = p['M_B']*p['x_BG']*p['y_BG']*u['omega_z']-p['M_B']*p['x_BG']*(p['y_BG']*u['omega_z']-2*p['z_BG']*qd[0])-p['M_B']*p['y_BG']*p['z_BG']*u['omega_x']+p['M_B']*p['z_BG']*(-2*p['x_BG']*qd[0]+p['y_BG']*u['omega_x'])
    return CC

def K_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear stiffness matrix 
    q:  degrees of freedom, array-like: ['phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['M_B', 'g', 'x_BG', 'z_BG']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    KK = np.zeros((1,1))
    KK[0,0] = p['M_B']*p['g']*p['x_BG']*sin(q[0])-p['M_B']*p['g']*p['z_BG']*cos(q[0])
    return KK

def B_lin(q=None,qd=None,p=None,u=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['phi_y(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: ['J_xx_B', 'J_zz_B', 'M_B', 'x_BG', 'y_BG', 'z_BG']
    u:  inputs at operating point, dictionary with keys: ['omega_x', 'omega_z']
           where each values is a constant!
    The columns of B correspond to:   [Derivative(omega_x(t), t), Derivative(omega_z(t), t), omega_x(t), omega_z(t)]\\ 
    """
    BB = np.zeros((1,4))
    BB[0,0] = 0
    BB[0,1] = 0
    BB[0,2] = 2*p['M_B']*p['x_BG']*p['z_BG']*u['omega_x']+p['M_B']*p['x_BG']*(p['x_BG']*u['omega_z']-2*p['z_BG']*u['omega_x'])+p['M_B']*p['y_BG']*p['z_BG']*qd[0]-p['M_B']*p['z_BG']*(p['y_BG']*qd[0]+p['z_BG']*u['omega_z'])+u['omega_z']*(-p['J_xx_B']+p['M_B']*(p['y_BG']**2+p['z_BG']**2))+u['omega_z']*(p['J_zz_B']-p['M_B']*(p['x_BG']**2+p['y_BG']**2))
    BB[0,3] = -p['M_B']*p['x_BG']*p['y_BG']*qd[0]-2*p['M_B']*p['x_BG']*p['z_BG']*u['omega_z']+p['M_B']*p['x_BG']*(p['x_BG']*u['omega_x']+p['y_BG']*qd[0])-p['M_B']*p['z_BG']*(-2*p['x_BG']*u['omega_z']+p['z_BG']*u['omega_x'])-u['omega_x']*(p['J_xx_B']-p['M_B']*(p['y_BG']**2+p['z_BG']**2))+u['omega_x']*(p['J_zz_B']-p['M_B']*(p['x_BG']**2+p['y_BG']**2))
    return BB

