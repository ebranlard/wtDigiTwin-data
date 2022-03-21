"""
Equations of motion
model name: F100010T0RNA_fnd_noLoads
"""
import numpy as np
from numpy import cos, sin, pi, sqrt
def info():
    """ Return information about current model present in this package """
    I=dict()
    I['name']='F100010T0RNA_fnd_noLoads'
    I['nq']=2
    I['nu']=0
    I['sq']=['x','phi_y']
    I['su']=[]
    return I

def forcing(t,q=None,qd=None,p=None,u=None,z=None):
    """ Non linear mass forcing 
    q:  degrees of freedom, array-like: ['x(t)', 'phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['L_T', 'M_F', 'M_RNA', 'M_T', 'g', 'x_RNAG', 'z_FG', 'z_RNAG', 'z_TG']
    u:  inputs, dictionary with keys: []
           where each values is a function of time
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    FF = np.zeros((2,1))
    FF[0,0] = p['L_T']*p['M_RNA']*qd[1]**2*sin(q[1])+p['M_F']*p['z_FG']*qd[1]**2*sin(q[1])+p['M_RNA']*p['x_RNAG']*qd[1]**2*cos(q[1])+p['M_RNA']*p['z_RNAG']*qd[1]**2*sin(q[1])+p['M_T']*p['z_TG']*qd[1]**2*sin(q[1])
    FF[1,0] = p['L_T']*p['M_RNA']*p['g']*sin(q[1])+p['M_F']*p['g']*p['z_FG']*sin(q[1])+p['M_RNA']*p['g']*p['x_RNAG']*cos(q[1])+p['M_RNA']*p['g']*p['z_RNAG']*sin(q[1])+p['M_T']*p['g']*p['z_TG']*sin(q[1])
    return FF

def mass_matrix(q=None,p=None,z=None):
    """ Non linear mass matrix 
     q:  degrees of freedom, array-like: ['x(t)', 'phi_y(t)']
     p:  parameters, dictionary with keys: ['J_yy_F', 'J_yy_RNA', 'J_yy_T', 'L_T', 'M_F', 'M_RNA', 'M_T', 'x_RNAG', 'z_FG', 'z_RNAG', 'z_TG']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((2,2))
    MM[0,0] = p['M_F']+p['M_RNA']+p['M_T']
    MM[0,1] = p['L_T']*p['M_RNA']*cos(q[1])+p['M_F']*p['z_FG']*cos(q[1])-p['M_RNA']*p['x_RNAG']*sin(q[1])+p['M_RNA']*p['z_RNAG']*cos(q[1])+p['M_T']*p['z_TG']*cos(q[1])
    MM[1,0] = p['L_T']*p['M_RNA']*cos(q[1])+p['M_F']*p['z_FG']*cos(q[1])-p['M_RNA']*p['x_RNAG']*sin(q[1])+p['M_RNA']*p['z_RNAG']*cos(q[1])+p['M_T']*p['z_TG']*cos(q[1])
    MM[1,1] = p['J_yy_F']+p['J_yy_RNA']+p['J_yy_T']+p['L_T']**2*p['M_RNA']+2*p['L_T']*p['M_RNA']*p['z_RNAG']+p['M_F']*p['z_FG']**2+p['M_RNA']*p['x_RNAG']**2+p['M_RNA']*p['z_RNAG']**2+p['M_T']*p['z_TG']**2
    return MM

def M_lin(q=None,p=None,z=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'phi_y(t)']
    p:  parameters, dictionary with keys: ['J_yy_F', 'J_yy_RNA', 'J_yy_T', 'L_T', 'M_F', 'M_RNA', 'M_T', 'x_RNAG', 'z_FG', 'z_RNAG', 'z_TG']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((2,2))
    MM[0,0] = p['M_F']+p['M_RNA']+p['M_T']
    MM[0,1] = p['M_F']*p['z_FG']*cos(q[1])+p['M_RNA']*(p['L_T']*cos(q[1])-p['x_RNAG']*sin(q[1])+p['z_RNAG']*cos(q[1]))+p['M_T']*p['z_TG']*cos(q[1])
    MM[1,0] = p['M_F']*p['z_FG']*cos(q[1])+p['M_RNA']*(p['L_T']*cos(q[1])-p['x_RNAG']*sin(q[1])+p['z_RNAG']*cos(q[1]))+p['M_T']*p['z_TG']*cos(q[1])
    MM[1,1] = p['J_yy_F']+p['J_yy_RNA']+p['J_yy_T']+p['M_F']*p['z_FG']**2+p['M_RNA']*(p['L_T']**2+2*p['L_T']*p['z_RNAG']+p['x_RNAG']**2+p['z_RNAG']**2)+p['M_T']*p['z_TG']**2
    return MM

def C_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear damping matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'phi_y(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: ['L_T', 'M_F', 'M_RNA', 'M_T', 'x_RNAG', 'z_FG', 'z_RNAG', 'z_TG']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    CC = np.zeros((2,2))
    CC[0,0] = 0
    CC[0,1] = -2*p['L_T']*p['M_RNA']*qd[1]*sin(q[1])-2*p['M_F']*p['z_FG']*qd[1]*sin(q[1])-2*p['M_RNA']*p['x_RNAG']*qd[1]*cos(q[1])-2*p['M_RNA']*p['z_RNAG']*qd[1]*sin(q[1])-2*p['M_T']*p['z_TG']*qd[1]*sin(q[1])
    CC[1,0] = 0
    CC[1,1] = 0
    return CC

def K_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear stiffness matrix 
    q:  degrees of freedom, array-like: ['x(t)', 'phi_y(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['L_T', 'M_F', 'M_RNA', 'M_T', 'g', 'x_RNAG', 'z_FG', 'z_RNAG', 'z_TG']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    KK = np.zeros((2,2))
    KK[0,0] = 0
    KK[0,1] = -p['L_T']*p['M_RNA']*qd[1]**2*cos(q[1])-p['M_F']*p['z_FG']*qd[1]**2*cos(q[1])+p['M_RNA']*p['x_RNAG']*qd[1]**2*sin(q[1])-p['M_RNA']*p['z_RNAG']*qd[1]**2*cos(q[1])-p['M_T']*p['z_TG']*qd[1]**2*cos(q[1])
    KK[1,0] = 0
    KK[1,1] = -p['L_T']*p['M_RNA']*p['g']*cos(q[1])-p['M_F']*p['g']*p['z_FG']*cos(q[1])+p['M_RNA']*p['g']*p['x_RNAG']*sin(q[1])-p['M_RNA']*p['g']*p['z_RNAG']*cos(q[1])-p['M_T']*p['g']*p['z_TG']*cos(q[1])
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

