"""
Equations of motion
model name: F100000T0RNA_fnd_noLoads
"""
import numpy as np
from numpy import cos, sin, pi, sqrt
def info():
    """ Return information about current model present in this package """
    I=dict()
    I['name']='F100000T0RNA_fnd_noLoads'
    I['nq']=1
    I['nu']=0
    I['sq']=['x']
    I['su']=[]
    return I

def forcing(t,q=None,qd=None,p=None,u=None,z=None):
    """ Non linear mass forcing 
    q:  degrees of freedom, array-like: ['x(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: []
    u:  inputs, dictionary with keys: []
           where each values is a function of time
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    FF = np.zeros((1,1))
    FF[0,0] = 0
    return FF

def mass_matrix(q=None,p=None,z=None):
    """ Non linear mass matrix 
     q:  degrees of freedom, array-like: ['x(t)']
     p:  parameters, dictionary with keys: ['M_F', 'M_RNA', 'M_T']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((1,1))
    MM[0,0] = p['M_F']+p['M_RNA']+p['M_T']
    return MM

def M_lin(q=None,p=None,z=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)']
    p:  parameters, dictionary with keys: ['M_F', 'M_RNA', 'M_T']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((1,1))
    MM[0,0] = p['M_F']+p['M_RNA']+p['M_T']
    return MM

def C_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear damping matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)']
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
    q:  degrees of freedom, array-like: ['x(t)']
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
    q:  degrees of freedom at operating point, array-like: ['x(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: []
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    The columns of B correspond to:   []\\ 
    """
    BB = np.zeros((0,0))
    return BB

