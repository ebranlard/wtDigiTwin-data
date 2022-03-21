"""
Equations of motion
model name: F100010T1N0S1_fnd
"""
import numpy as np
from numpy import cos, sin, pi, sqrt
def info():
    """ Return information about current model present in this package """
    I=dict()
    I['name']='F100010T1N0S1_fnd'
    I['nq']=4
    I['nu']=2
    I['sq']=['x','phi_y','q_T1','psi']
    I['su']=['F_hx','F_hz']
    return I

def forcing(t,q=None,qd=None,p=None,u=None,z=None):
    """ Non linear mass forcing 
    q:  degrees of freedom, array-like: ['x(t)', 'phi_y(t)', 'q_T1(t)', 'psi(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['DD_T', 'KK_T', 'KM_00', 'KM_04', 'KM_44', 'L_T', 'MM_T', 'M_F', 'M_N', 'M_R', 'g', 'v_yT1c', 'x_NG', 'x_NR', 'z_FG', 'z_NG', 'z_NR', 'z_T0']
    u:  inputs, dictionary with keys: ['F_hx', 'F_hz']
           where each values is a function of time
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    FF = np.zeros((4,1))
    FF[0,0] = 2*p['MM_T'][0,6]*qd[1]*qd[2]*sin(q[1])-p['MM_T'][1,3]*qd[1]**2*sin(q[1])-p['KM_00']*q[0]-p['KM_04']*q[1]+p['L_T']*p['M_N']*qd[1]**2*sin(q[1])+p['L_T']*p['M_R']*qd[1]**2*sin(q[1])+p['M_F']*p['z_FG']*qd[1]**2*sin(q[1])-p['M_N']*p['v_yT1c']**3*p['x_NG']*q[2]*qd[2]**2*sin(q[1])+p['M_N']*p['v_yT1c']**3*p['z_NG']*q[2]*qd[2]**2*cos(q[1])-2*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[2]*qd[1]*qd[2]*sin(q[1])+p['M_N']*p['v_yT1c']**2*p['x_NG']*qd[2]**2*cos(q[1])+2*p['M_N']*p['v_yT1c']**2*p['z_NG']*q[2]*qd[1]*qd[2]*cos(q[1])+p['M_N']*p['v_yT1c']**2*p['z_NG']*qd[2]**2*sin(q[1])-p['M_N']*p['v_yT1c']*p['x_NG']*q[2]*qd[1]**2*sin(q[1])+2*p['M_N']*p['v_yT1c']*p['x_NG']*qd[1]*qd[2]*cos(q[1])+p['M_N']*p['v_yT1c']*p['z_NG']*q[2]*qd[1]**2*cos(q[1])+2*p['M_N']*p['v_yT1c']*p['z_NG']*qd[1]*qd[2]*sin(q[1])+p['M_N']*p['x_NG']*qd[1]**2*cos(q[1])+p['M_N']*p['z_NG']*qd[1]**2*sin(q[1])+p['M_N']*q[2]*qd[1]**2*cos(q[1])+2*p['M_N']*qd[1]*qd[2]*sin(q[1])-p['M_R']*p['v_yT1c']**3*p['x_NR']*q[2]*qd[2]**2*sin(q[1])+p['M_R']*p['v_yT1c']**3*p['z_NR']*q[2]*qd[2]**2*cos(q[1])-2*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[2]*qd[1]*qd[2]*sin(q[1])+p['M_R']*p['v_yT1c']**2*p['x_NR']*qd[2]**2*cos(q[1])+2*p['M_R']*p['v_yT1c']**2*p['z_NR']*q[2]*qd[1]*qd[2]*cos(q[1])+p['M_R']*p['v_yT1c']**2*p['z_NR']*qd[2]**2*sin(q[1])-p['M_R']*p['v_yT1c']*p['x_NR']*q[2]*qd[1]**2*sin(q[1])+2*p['M_R']*p['v_yT1c']*p['x_NR']*qd[1]*qd[2]*cos(q[1])+p['M_R']*p['v_yT1c']*p['z_NR']*q[2]*qd[1]**2*cos(q[1])+2*p['M_R']*p['v_yT1c']*p['z_NR']*qd[1]*qd[2]*sin(q[1])+p['M_R']*p['x_NR']*qd[1]**2*cos(q[1])+p['M_R']*p['z_NR']*qd[1]**2*sin(q[1])+p['M_R']*q[2]*qd[1]**2*cos(q[1])+2*p['M_R']*qd[1]*qd[2]*sin(q[1])+u['F_hx'](t,q,qd)
    FF[1,0] = -p['MM_T'][1,3]*p['g']*sin(q[1])-p['KM_04']*q[0]-p['KM_44']*q[1]+p['L_T']*p['M_N']*p['g']*sin(q[1])+p['L_T']*p['M_N']*p['v_yT1c']**3*p['z_NG']*q[2]*qd[2]**2+p['L_T']*p['M_N']*p['v_yT1c']**2*p['x_NG']*qd[2]**2+2*p['L_T']*p['M_N']*p['v_yT1c']**2*p['z_NG']*q[2]*qd[1]*qd[2]+2*p['L_T']*p['M_N']*p['v_yT1c']*p['x_NG']*qd[1]*qd[2]+p['L_T']*p['M_R']*p['g']*sin(q[1])+p['L_T']*p['M_R']*p['v_yT1c']**3*p['z_NR']*q[2]*qd[2]**2+p['L_T']*p['M_R']*p['v_yT1c']**2*p['x_NR']*qd[2]**2+2*p['L_T']*p['M_R']*p['v_yT1c']**2*p['z_NR']*q[2]*qd[1]*qd[2]+2*p['L_T']*p['M_R']*p['v_yT1c']*p['x_NR']*qd[1]*qd[2]+p['M_F']*p['g']*p['z_FG']*sin(q[1])-p['M_N']*p['g']*p['v_yT1c']*p['x_NG']*q[2]*sin(q[1])+p['M_N']*p['g']*p['v_yT1c']*p['z_NG']*q[2]*cos(q[1])+p['M_N']*p['g']*p['x_NG']*cos(q[1])+p['M_N']*p['g']*p['z_NG']*sin(q[1])+p['M_N']*p['g']*q[2]*cos(q[1])+p['M_N']*p['v_yT1c']**3*p['x_NG']*q[2]**2*qd[2]**2+2*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[2]**2*qd[1]*qd[2]-p['M_N']*p['v_yT1c']**2*p['z_NG']*q[2]*qd[2]**2-4*p['M_N']*p['v_yT1c']*p['z_NG']*q[2]*qd[1]*qd[2]-2*p['M_N']*p['x_NG']*qd[1]*qd[2]-2*p['M_N']*q[2]*qd[1]*qd[2]-p['M_R']*p['g']*p['v_yT1c']*p['x_NR']*q[2]*sin(q[1])+p['M_R']*p['g']*p['v_yT1c']*p['z_NR']*q[2]*cos(q[1])+p['M_R']*p['g']*p['x_NR']*cos(q[1])+p['M_R']*p['g']*p['z_NR']*sin(q[1])+p['M_R']*p['g']*q[2]*cos(q[1])+p['M_R']*p['v_yT1c']**3*p['x_NR']*q[2]**2*qd[2]**2+2*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[2]**2*qd[1]*qd[2]-p['M_R']*p['v_yT1c']**2*p['z_NR']*q[2]*qd[2]**2-4*p['M_R']*p['v_yT1c']*p['z_NR']*q[2]*qd[1]*qd[2]-2*p['M_R']*p['x_NR']*qd[1]*qd[2]-2*p['M_R']*q[2]*qd[1]*qd[2]+p['z_T0']*u['F_hx'](t,q,qd)*cos(q[1])-p['z_T0']*u['F_hz'](t,q,qd)*sin(q[1])
    FF[2,0] = p['MM_T'][0,6]*p['g']*sin(q[1])-p['DD_T'][6,6]*qd[2]-p['KK_T'][6,6]*q[2]-p['L_T']*p['M_N']*p['v_yT1c']**2*p['z_NG']*q[2]*qd[1]**2-p['L_T']*p['M_N']*p['v_yT1c']*p['x_NG']*qd[1]**2-p['L_T']*p['M_R']*p['v_yT1c']**2*p['z_NR']*q[2]*qd[1]**2-p['L_T']*p['M_R']*p['v_yT1c']*p['x_NR']*qd[1]**2-p['M_N']*p['g']*p['v_yT1c']**2*p['x_NG']*q[2]*sin(q[1])+p['M_N']*p['g']*p['v_yT1c']**2*p['z_NG']*q[2]*cos(q[1])+p['M_N']*p['g']*p['v_yT1c']*p['x_NG']*cos(q[1])+p['M_N']*p['g']*p['v_yT1c']*p['z_NG']*sin(q[1])+p['M_N']*p['g']*sin(q[1])+p['M_N']*p['v_yT1c']**3*p['z_NG']*q[2]*qd[2]**2-p['M_N']*p['v_yT1c']**2*p['x_NG']*q[2]**2*qd[1]**2+p['M_N']*p['v_yT1c']**2*p['x_NG']*qd[2]**2+2*p['M_N']*p['v_yT1c']*p['z_NG']*q[2]*qd[1]**2+p['M_N']*p['x_NG']*qd[1]**2+p['M_N']*q[2]*qd[1]**2-p['M_R']*p['g']*p['v_yT1c']**2*p['x_NR']*q[2]*sin(q[1])+p['M_R']*p['g']*p['v_yT1c']**2*p['z_NR']*q[2]*cos(q[1])+p['M_R']*p['g']*p['v_yT1c']*p['x_NR']*cos(q[1])+p['M_R']*p['g']*p['v_yT1c']*p['z_NR']*sin(q[1])+p['M_R']*p['g']*sin(q[1])+p['M_R']*p['v_yT1c']**3*p['z_NR']*q[2]*qd[2]**2-p['M_R']*p['v_yT1c']**2*p['x_NR']*q[2]**2*qd[1]**2+p['M_R']*p['v_yT1c']**2*p['x_NR']*qd[2]**2+2*p['M_R']*p['v_yT1c']*p['z_NR']*q[2]*qd[1]**2+p['M_R']*p['x_NR']*qd[1]**2+p['M_R']*q[2]*qd[1]**2
    FF[3,0] = 0
    return FF

def mass_matrix(q=None,p=None,z=None):
    """ Non linear mass matrix 
     q:  degrees of freedom, array-like: ['x(t)', 'phi_y(t)', 'q_T1(t)', 'psi(t)']
     p:  parameters, dictionary with keys: ['JO_R', 'J_yy_F', 'J_yy_N', 'Jxx_R', 'L_T', 'MM_T', 'M_F', 'M_N', 'M_R', 'v_yT1c', 'x_NG', 'x_NR', 'z_FG', 'z_NG', 'z_NR']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((4,4))
    MM[0,0] = p['MM_T'][0,0]*sin(q[1])**2+p['MM_T'][0,0]*cos(q[1])**2+p['M_F']+p['M_N']+p['M_R']
    MM[0,1] = -p['MM_T'][1,3]*cos(q[1])+p['L_T']*p['M_N']*cos(q[1])+p['L_T']*p['M_R']*cos(q[1])+p['M_F']*p['z_FG']*cos(q[1])-p['M_N']*p['v_yT1c']*p['x_NG']*q[2]*cos(q[1])-p['M_N']*p['v_yT1c']*p['z_NG']*q[2]*sin(q[1])-p['M_N']*p['x_NG']*sin(q[1])+p['M_N']*p['z_NG']*cos(q[1])-p['M_N']*q[2]*sin(q[1])-p['M_R']*p['v_yT1c']*p['x_NR']*q[2]*cos(q[1])-p['M_R']*p['v_yT1c']*p['z_NR']*q[2]*sin(q[1])-p['M_R']*p['x_NR']*sin(q[1])+p['M_R']*p['z_NR']*cos(q[1])-p['M_R']*q[2]*sin(q[1])
    MM[0,2] = p['MM_T'][0,6]*cos(q[1])-p['M_N']*p['v_yT1c']**2*p['x_NG']*q[2]*cos(q[1])-p['M_N']*p['v_yT1c']**2*p['z_NG']*q[2]*sin(q[1])-p['M_N']*p['v_yT1c']*p['x_NG']*sin(q[1])+p['M_N']*p['v_yT1c']*p['z_NG']*cos(q[1])+p['M_N']*cos(q[1])-p['M_R']*p['v_yT1c']**2*p['x_NR']*q[2]*cos(q[1])-p['M_R']*p['v_yT1c']**2*p['z_NR']*q[2]*sin(q[1])-p['M_R']*p['v_yT1c']*p['x_NR']*sin(q[1])+p['M_R']*p['v_yT1c']*p['z_NR']*cos(q[1])+p['M_R']*cos(q[1])
    MM[0,3] = 0
    MM[1,0] = -p['MM_T'][1,3]*cos(q[1])+p['L_T']*p['M_N']*cos(q[1])+p['L_T']*p['M_R']*cos(q[1])+p['M_F']*p['z_FG']*cos(q[1])-p['M_N']*p['v_yT1c']*p['x_NG']*q[2]*cos(q[1])-p['M_N']*p['v_yT1c']*p['z_NG']*q[2]*sin(q[1])-p['M_N']*p['x_NG']*sin(q[1])+p['M_N']*p['z_NG']*cos(q[1])-p['M_N']*q[2]*sin(q[1])-p['M_R']*p['v_yT1c']*p['x_NR']*q[2]*cos(q[1])-p['M_R']*p['v_yT1c']*p['z_NR']*q[2]*sin(q[1])-p['M_R']*p['x_NR']*sin(q[1])+p['M_R']*p['z_NR']*cos(q[1])-p['M_R']*q[2]*sin(q[1])
    MM[1,1] = p['MM_T'][4,4]+p['JO_R']*sin(q[3])**2+p['JO_R']*cos(q[3])**2+p['J_yy_F']+p['J_yy_N']+p['L_T']**2*p['M_N']+p['L_T']**2*p['M_R']-2*p['L_T']*p['M_N']*p['v_yT1c']*p['x_NG']*q[2]+2*p['L_T']*p['M_N']*p['z_NG']-2*p['L_T']*p['M_R']*p['v_yT1c']*p['x_NR']*q[2]+2*p['L_T']*p['M_R']*p['z_NR']+p['M_F']*p['z_FG']**2+2*p['M_N']*p['v_yT1c']*p['z_NG']*q[2]**2+p['M_N']*p['x_NG']**2+2*p['M_N']*p['x_NG']*q[2]+p['M_N']*p['z_NG']**2+p['M_N']*q[2]**2+2*p['M_R']*p['v_yT1c']*p['z_NR']*q[2]**2+p['M_R']*p['x_NR']**2+2*p['M_R']*p['x_NR']*q[2]+p['M_R']*p['z_NR']**2+p['M_R']*q[2]**2
    MM[1,2] = p['MM_T'][6,4]+p['JO_R']*p['v_yT1c']*sin(q[3])**2+p['JO_R']*p['v_yT1c']*cos(q[3])**2+p['J_yy_N']*p['v_yT1c']-p['L_T']*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[2]+p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']+p['L_T']*p['M_N']-p['L_T']*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[2]+p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']+p['L_T']*p['M_R']+p['M_N']*p['v_yT1c']**2*p['z_NG']*q[2]**2+p['M_N']*p['v_yT1c']*p['x_NG']**2+p['M_N']*p['v_yT1c']*p['z_NG']**2+p['M_N']*p['z_NG']+p['M_R']*p['v_yT1c']**2*p['z_NR']*q[2]**2+p['M_R']*p['v_yT1c']*p['x_NR']**2+p['M_R']*p['v_yT1c']*p['z_NR']**2+p['M_R']*p['z_NR']
    MM[1,3] = 0
    MM[2,0] = p['MM_T'][0,6]*cos(q[1])-p['M_N']*p['v_yT1c']**2*p['x_NG']*q[2]*cos(q[1])-p['M_N']*p['v_yT1c']**2*p['z_NG']*q[2]*sin(q[1])-p['M_N']*p['v_yT1c']*p['x_NG']*sin(q[1])+p['M_N']*p['v_yT1c']*p['z_NG']*cos(q[1])+p['M_N']*cos(q[1])-p['M_R']*p['v_yT1c']**2*p['x_NR']*q[2]*cos(q[1])-p['M_R']*p['v_yT1c']**2*p['z_NR']*q[2]*sin(q[1])-p['M_R']*p['v_yT1c']*p['x_NR']*sin(q[1])+p['M_R']*p['v_yT1c']*p['z_NR']*cos(q[1])+p['M_R']*cos(q[1])
    MM[2,1] = p['MM_T'][6,4]+p['JO_R']*p['v_yT1c']*sin(q[3])**2+p['JO_R']*p['v_yT1c']*cos(q[3])**2+p['J_yy_N']*p['v_yT1c']-p['L_T']*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[2]+p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']+p['L_T']*p['M_N']-p['L_T']*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[2]+p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']+p['L_T']*p['M_R']+p['M_N']*p['v_yT1c']**2*p['z_NG']*q[2]**2+p['M_N']*p['v_yT1c']*p['x_NG']**2+p['M_N']*p['v_yT1c']*p['z_NG']**2+p['M_N']*p['z_NG']+p['M_R']*p['v_yT1c']**2*p['z_NR']*q[2]**2+p['M_R']*p['v_yT1c']*p['x_NR']**2+p['M_R']*p['v_yT1c']*p['z_NR']**2+p['M_R']*p['z_NR']
    MM[2,2] = p['MM_T'][6,6]+p['JO_R']*p['v_yT1c']**2*sin(q[3])**2+p['JO_R']*p['v_yT1c']**2*cos(q[3])**2+p['J_yy_N']*p['v_yT1c']**2+p['M_N']*p['v_yT1c']**2*p['x_NG']**2-2*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[2]+p['M_N']*p['v_yT1c']**2*p['z_NG']**2+2*p['M_N']*p['v_yT1c']*p['z_NG']+p['M_N']+p['M_R']*p['v_yT1c']**2*p['x_NR']**2-2*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[2]+p['M_R']*p['v_yT1c']**2*p['z_NR']**2+2*p['M_R']*p['v_yT1c']*p['z_NR']+p['M_R']
    MM[2,3] = 0
    MM[3,0] = 0
    MM[3,1] = 0
    MM[3,2] = 0
    MM[3,3] = p['Jxx_R']
    return MM

def M_lin(q=None,p=None,z=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'phi_y(t)', 'q_T1(t)', 'psi(t)']
    p:  parameters, dictionary with keys: ['JO_R', 'J_yy_F', 'J_yy_N', 'Jxx_R', 'L_T', 'MM_T', 'M_F', 'M_N', 'M_R', 'v_yT1c', 'x_NG', 'x_NR', 'z_FG', 'z_NG', 'z_NR']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((4,4))
    MM[0,0] = p['MM_T'][0,0]*sin(q[1])**2+p['MM_T'][0,0]*cos(q[1])**2+p['M_F']+p['M_N']+p['M_R']
    MM[0,1] = -p['MM_T'][1,3]*cos(q[1])+p['M_F']*p['z_FG']*cos(q[1])+p['M_N']*(p['L_T']*cos(q[1])-p['x_NG']*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))+p['z_NG']*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-q[2]*sin(q[1]))+p['M_R']*(p['L_T']*cos(q[1])-p['x_NR']*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))+p['z_NR']*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-q[2]*sin(q[1]))
    MM[0,2] = p['MM_T'][0,6]*cos(q[1])+p['M_N']*(-p['v_yT1c']*p['x_NG']*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))+p['v_yT1c']*p['z_NG']*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))+cos(q[1]))+p['M_R']*(-p['v_yT1c']*p['x_NR']*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))+p['v_yT1c']*p['z_NR']*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))+cos(q[1]))
    MM[0,3] = 0
    MM[1,0] = -p['MM_T'][1,3]*cos(q[1])+p['M_F']*p['z_FG']*cos(q[1])+p['M_N']*(p['L_T']*cos(q[1])-p['x_NG']*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))+p['z_NG']*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-q[2]*sin(q[1]))+p['M_R']*(p['L_T']*cos(q[1])-p['x_NR']*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))+p['z_NR']*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-q[2]*sin(q[1]))
    MM[1,1] = p['MM_T'][4,4]+p['JO_R']*sin(q[3])**2+p['JO_R']*cos(q[3])**2+p['J_yy_F']+p['J_yy_N']+p['M_F']*p['z_FG']**2+p['M_N']*(p['L_T']**2+p['L_T']*(-p['v_yT1c']*p['x_NG']*q[2]+p['z_NG'])+p['x_NG']**2-p['x_NG']*(p['L_T']*p['v_yT1c']*q[2]-q[2])+p['z_NG']**2+p['z_NG']*(p['L_T']+p['v_yT1c']*q[2]**2)+q[2]**2-q[2]*(-p['v_yT1c']*p['z_NG']*q[2]-p['x_NG']))+p['M_R']*(p['L_T']**2+p['L_T']*(-p['v_yT1c']*p['x_NR']*q[2]+p['z_NR'])+p['x_NR']**2-p['x_NR']*(p['L_T']*p['v_yT1c']*q[2]-q[2])+p['z_NR']**2+p['z_NR']*(p['L_T']+p['v_yT1c']*q[2]**2)+q[2]**2-q[2]*(-p['v_yT1c']*p['z_NR']*q[2]-p['x_NR']))
    MM[1,2] = p['MM_T'][6,4]+p['JO_R']*p['v_yT1c']*sin(q[3])**2+p['JO_R']*p['v_yT1c']*cos(q[3])**2+p['J_yy_N']*p['v_yT1c']+p['M_N']*(p['L_T']*(-p['v_yT1c']**2*p['x_NG']*q[2]+p['v_yT1c']*p['z_NG'])+p['L_T']+p['v_yT1c']*p['x_NG']**2-p['v_yT1c']*p['x_NG']*q[2]+p['v_yT1c']*p['z_NG']**2+p['z_NG']-q[2]*(-p['v_yT1c']**2*p['z_NG']*q[2]-p['v_yT1c']*p['x_NG']))+p['M_R']*(p['L_T']*(-p['v_yT1c']**2*p['x_NR']*q[2]+p['v_yT1c']*p['z_NR'])+p['L_T']+p['v_yT1c']*p['x_NR']**2-p['v_yT1c']*p['x_NR']*q[2]+p['v_yT1c']*p['z_NR']**2+p['z_NR']-q[2]*(-p['v_yT1c']**2*p['z_NR']*q[2]-p['v_yT1c']*p['x_NR']))
    MM[1,3] = 0
    MM[2,0] = p['MM_T'][0,6]*cos(q[1])+p['M_N']*(-p['v_yT1c']*p['x_NG']*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))+p['v_yT1c']*p['z_NG']*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))+cos(q[1]))+p['M_R']*(-p['v_yT1c']*p['x_NR']*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))+p['v_yT1c']*p['z_NR']*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))+cos(q[1]))
    MM[2,1] = p['MM_T'][6,4]+p['JO_R']*p['v_yT1c']*sin(q[3])**2+p['JO_R']*p['v_yT1c']*cos(q[3])**2+p['J_yy_N']*p['v_yT1c']+p['M_N']*(p['L_T']+p['v_yT1c']*p['x_NG']**2-p['v_yT1c']*p['x_NG']*q[2]-p['v_yT1c']*p['x_NG']*(p['L_T']*p['v_yT1c']*q[2]-q[2])+p['v_yT1c']*p['z_NG']**2+p['v_yT1c']*p['z_NG']*(p['L_T']+p['v_yT1c']*q[2]**2)+p['z_NG'])+p['M_R']*(p['L_T']+p['v_yT1c']*p['x_NR']**2-p['v_yT1c']*p['x_NR']*q[2]-p['v_yT1c']*p['x_NR']*(p['L_T']*p['v_yT1c']*q[2]-q[2])+p['v_yT1c']*p['z_NR']**2+p['v_yT1c']*p['z_NR']*(p['L_T']+p['v_yT1c']*q[2]**2)+p['z_NR'])
    MM[2,2] = p['MM_T'][6,6]+p['JO_R']*p['v_yT1c']**2*sin(q[3])**2+p['JO_R']*p['v_yT1c']**2*cos(q[3])**2+p['J_yy_N']*p['v_yT1c']**2+p['M_N']*(p['v_yT1c']**2*p['x_NG']**2-2*p['v_yT1c']**2*p['x_NG']*q[2]+p['v_yT1c']**2*p['z_NG']**2+2*p['v_yT1c']*p['z_NG']+1)+p['M_R']*(p['v_yT1c']**2*p['x_NR']**2-2*p['v_yT1c']**2*p['x_NR']*q[2]+p['v_yT1c']**2*p['z_NR']**2+2*p['v_yT1c']*p['z_NR']+1)
    MM[2,3] = 0
    MM[3,0] = 0
    MM[3,1] = 0
    MM[3,2] = 0
    MM[3,3] = p['Jxx_R']
    return MM

def C_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear damping matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'phi_y(t)', 'q_T1(t)', 'psi(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: ['DD_T', 'L_T', 'MM_T', 'M_F', 'M_N', 'M_R', 'v_yT1c', 'x_NG', 'x_NR', 'z_FG', 'z_NG', 'z_NR']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    CC = np.zeros((4,4))
    CC[0,0] = 0
    CC[0,1] = -2*p['M_F']*p['z_FG']*qd[1]*sin(q[1])-p['M_N']*p['x_NG']*(p['v_yT1c']*qd[2]+qd[1])*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-p['M_N']*p['z_NG']*(p['v_yT1c']*qd[2]+qd[1])*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))-2*p['M_N']*q[2]*qd[1]*cos(q[1])+p['M_N']*(-2*p['L_T']*qd[1]-2*qd[2])*sin(q[1])+p['M_N']*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-p['M_N']*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))-p['M_R']*p['x_NR']*(2*p['v_yT1c']*qd[2]+2*qd[1])*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-p['M_R']*p['z_NR']*(2*p['v_yT1c']*qd[2]+2*qd[1])*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))-2*p['M_R']*q[2]*qd[1]*cos(q[1])+p['M_R']*(-2*p['L_T']*qd[1]-2*qd[2])*sin(q[1])-(2*p['MM_T'][0,6]*qd[2]-2*p['MM_T'][1,3]*qd[1])*sin(q[1])
    CC[0,2] = -2*p['MM_T'][0,6]*qd[1]*sin(q[1])-p['M_N']*p['v_yT1c']*p['x_NG']*(p['v_yT1c']*qd[2]+qd[1])*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-p['M_N']*p['v_yT1c']*p['z_NG']*(p['v_yT1c']*qd[2]+qd[1])*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))+p['M_N']*p['v_yT1c']*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-p['M_N']*p['v_yT1c']*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))-2*p['M_N']*qd[1]*sin(q[1])-2*p['M_R']*p['v_yT1c']*p['x_NR']*(p['v_yT1c']*qd[2]+qd[1])*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))-2*p['M_R']*p['v_yT1c']*p['z_NR']*(p['v_yT1c']*qd[2]+qd[1])*(p['v_yT1c']*q[2]*cos(q[1])+sin(q[1]))-2*p['M_R']*qd[1]*sin(q[1])
    CC[0,3] = 0
    CC[1,0] = 0
    CC[1,1] = -2*p['L_T']*p['M_N']*q[2]*qd[1]-2*p['L_T']*p['M_R']*q[2]*qd[1]-p['M_N']*p['x_NG']*(p['L_T']+p['v_yT1c']*q[2]**2)*(p['v_yT1c']*qd[2]+qd[1])+p['M_N']*p['x_NG']*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])-p['M_N']*p['z_NG']*(p['v_yT1c']*qd[2]+qd[1])*(p['L_T']*p['v_yT1c']*q[2]-q[2])+p['M_N']*p['z_NG']*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])-2*p['M_N']*q[2]*qd[1]*(-p['v_yT1c']*p['x_NG']*q[2]+p['z_NG'])-p['M_N']*q[2]*(-2*p['L_T']*qd[1]-2*qd[2])+p['M_N']*(p['L_T']+p['v_yT1c']*q[2]**2)*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])+p['M_N']*(-2*p['L_T']*qd[1]-2*qd[2])*(-p['v_yT1c']*p['z_NG']*q[2]-p['x_NG'])-p['M_N']*(p['L_T']*p['v_yT1c']*q[2]-q[2])*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])-p['M_R']*p['x_NR']*(p['L_T']+p['v_yT1c']*q[2]**2)*(2*p['v_yT1c']*qd[2]+2*qd[1])-p['M_R']*p['z_NR']*(2*p['v_yT1c']*qd[2]+2*qd[1])*(p['L_T']*p['v_yT1c']*q[2]-q[2])-2*p['M_R']*q[2]*qd[1]*(-p['v_yT1c']*p['x_NR']*q[2]+p['z_NR'])-p['M_R']*q[2]*(-2*p['L_T']*qd[1]-2*qd[2])+p['M_R']*(-2*p['L_T']*qd[1]-2*qd[2])*(-p['v_yT1c']*p['z_NR']*q[2]-p['x_NR'])
    CC[1,2] = -p['M_N']*p['v_yT1c']*p['x_NG']*(p['L_T']+p['v_yT1c']*q[2]**2)*(p['v_yT1c']*qd[2]+qd[1])+p['M_N']*p['v_yT1c']*p['x_NG']*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])-p['M_N']*p['v_yT1c']*p['z_NG']*(p['v_yT1c']*qd[2]+qd[1])*(p['L_T']*p['v_yT1c']*q[2]-q[2])+p['M_N']*p['v_yT1c']*p['z_NG']*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])+p['M_N']*p['v_yT1c']*(p['L_T']+p['v_yT1c']*q[2]**2)*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])-p['M_N']*p['v_yT1c']*(p['L_T']*p['v_yT1c']*q[2]-q[2])*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])+2*p['M_N']*q[2]*qd[1]-2*p['M_N']*qd[1]*(-p['v_yT1c']*p['z_NG']*q[2]-p['x_NG'])-2*p['M_R']*p['v_yT1c']*p['x_NR']*(p['L_T']+p['v_yT1c']*q[2]**2)*(p['v_yT1c']*qd[2]+qd[1])-2*p['M_R']*p['v_yT1c']*p['z_NR']*(p['v_yT1c']*qd[2]+qd[1])*(p['L_T']*p['v_yT1c']*q[2]-q[2])+2*p['M_R']*q[2]*qd[1]-2*p['M_R']*qd[1]*(-p['v_yT1c']*p['z_NR']*q[2]-p['x_NR'])
    CC[1,3] = 0
    CC[2,0] = 0
    CC[2,1] = p['M_N']*p['v_yT1c']*p['x_NG']*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])-p['M_N']*p['v_yT1c']*p['z_NG']*q[2]*(p['v_yT1c']*qd[2]+qd[1])+p['M_N']*p['v_yT1c']*p['z_NG']*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])-p['M_N']*p['v_yT1c']*q[2]*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])-p['M_N']*p['x_NG']*(p['v_yT1c']*qd[2]+qd[1])-2*p['M_N']*q[2]*qd[1]*(-p['v_yT1c']**2*p['x_NG']*q[2]+p['v_yT1c']*p['z_NG'])-2*p['M_N']*q[2]*qd[1]+p['M_N']*(-2*p['L_T']*qd[1]-2*qd[2])*(-p['v_yT1c']**2*p['z_NG']*q[2]-p['v_yT1c']*p['x_NG'])+p['M_N']*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])-p['M_R']*p['v_yT1c']*p['z_NR']*q[2]*(2*p['v_yT1c']*qd[2]+2*qd[1])-p['M_R']*p['x_NR']*(2*p['v_yT1c']*qd[2]+2*qd[1])-2*p['M_R']*q[2]*qd[1]*(-p['v_yT1c']**2*p['x_NR']*q[2]+p['v_yT1c']*p['z_NR'])-2*p['M_R']*q[2]*qd[1]+p['M_R']*(-2*p['L_T']*qd[1]-2*qd[2])*(-p['v_yT1c']**2*p['z_NR']*q[2]-p['v_yT1c']*p['x_NR'])
    CC[2,2] = p['DD_T'][6,6]+p['M_N']*p['v_yT1c']**2*p['x_NG']*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])-p['M_N']*p['v_yT1c']**2*p['z_NG']*q[2]*(p['v_yT1c']*qd[2]+qd[1])+p['M_N']*p['v_yT1c']**2*p['z_NG']*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])-p['M_N']*p['v_yT1c']**2*q[2]*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])-p['M_N']*p['v_yT1c']*p['x_NG']*(p['v_yT1c']*qd[2]+qd[1])+p['M_N']*p['v_yT1c']*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])-2*p['M_N']*qd[1]*(-p['v_yT1c']**2*p['z_NG']*q[2]-p['v_yT1c']*p['x_NG'])-2*p['M_R']*p['v_yT1c']**2*p['z_NR']*q[2]*(p['v_yT1c']*qd[2]+qd[1])-2*p['M_R']*p['v_yT1c']*p['x_NR']*(p['v_yT1c']*qd[2]+qd[1])-2*p['M_R']*qd[1]*(-p['v_yT1c']**2*p['z_NR']*q[2]-p['v_yT1c']*p['x_NR'])
    CC[2,3] = 0
    CC[3,0] = 0
    CC[3,1] = 0
    CC[3,2] = 0
    CC[3,3] = 0
    return CC

def K_lin(q=None,qd=None,p=None,u=None,z=None):
    """ Linear stiffness matrix 
    q:  degrees of freedom, array-like: ['x(t)', 'phi_y(t)', 'q_T1(t)', 'psi(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['KK_T', 'KM_00', 'KM_04', 'KM_44', 'L_T', 'MM_T', 'M_F', 'M_N', 'M_R', 'g', 'v_yT1c', 'x_NG', 'x_NR', 'z_FG', 'z_NG', 'z_NR', 'z_T0']
    u:  inputs at operating point, dictionary with keys: ['F_hx', 'F_hz']
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    KK = np.zeros((4,4))
    KK[0,0] = p['KM_00']
    KK[0,1] = -p['MM_T'][0,0]*p['g']*cos(q[1])**2+p['KM_04']-p['M_F']*p['z_FG']*qd[1]**2*cos(q[1])+p['M_N']*q[2]*qd[1]**2*sin(q[1])+p['M_N']*(p['v_yT1c']*qd[2]+qd[1])*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])*(-p['v_yT1c']*q[2]*cos(q[1])-sin(q[1]))-p['M_N']*(p['v_yT1c']*qd[2]+qd[1])*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))+p['M_N']*(-qd[1]*qd[2]-qd[1]*(p['L_T']*qd[1]+qd[2]))*cos(q[1])-p['M_R']*p['x_NR']*(p['v_yT1c']*qd[2]+qd[1])**2*(-p['v_yT1c']*q[2]*cos(q[1])-sin(q[1]))-p['M_R']*p['z_NR']*(p['v_yT1c']*qd[2]+qd[1])**2*(-p['v_yT1c']*q[2]*sin(q[1])+cos(q[1]))+p['M_R']*q[2]*qd[1]**2*sin(q[1])+p['M_R']*(-qd[1]*qd[2]-qd[1]*(p['L_T']*qd[1]+qd[2]))*cos(q[1])-(2*p['MM_T'][0,6]*qd[1]*qd[2]-p['MM_T'][0,0]*p['g']*cos(q[1])-p['MM_T'][1,3]*qd[1]**2)*cos(q[1])
    KK[0,2] = -p['M_N']*p['v_yT1c']*(p['v_yT1c']*qd[2]+qd[1])*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])*sin(q[1])-p['M_N']*p['v_yT1c']*(p['v_yT1c']*qd[2]+qd[1])*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])*cos(q[1])-p['M_N']*qd[1]**2*cos(q[1])+p['M_R']*p['v_yT1c']*p['x_NR']*(p['v_yT1c']*qd[2]+qd[1])**2*sin(q[1])-p['M_R']*p['v_yT1c']*p['z_NR']*(p['v_yT1c']*qd[2]+qd[1])**2*cos(q[1])-p['M_R']*qd[1]**2*cos(q[1])
    KK[0,3] = 0
    KK[1,0] = p['KM_04']
    KK[1,1] = p['MM_T'][1,3]*p['g']*cos(q[1])+p['KM_44']-p['L_T']*p['M_N']*p['g']*cos(q[1])-p['L_T']*p['M_R']*p['g']*cos(q[1])-p['M_F']*p['g']*p['z_FG']*cos(q[1])-p['M_N']*p['g']*p['x_NG']*(-p['v_yT1c']*q[2]*cos(q[1])-sin(q[1]))+p['M_N']*p['g']*p['z_NG']*(p['v_yT1c']*q[2]*sin(q[1])-cos(q[1]))+p['M_N']*p['g']*q[2]*sin(q[1])-p['M_R']*p['g']*p['x_NR']*(-p['v_yT1c']*q[2]*cos(q[1])-sin(q[1]))+p['M_R']*p['g']*p['z_NR']*(p['v_yT1c']*q[2]*sin(q[1])-cos(q[1]))+p['M_R']*p['g']*q[2]*sin(q[1])-p['z_T0']*(-u['F_hx']*sin(q[1])-u['F_hz']*cos(q[1]))
    KK[1,2] = -p['L_T']*p['M_N']*qd[1]**2-p['L_T']*p['M_R']*qd[1]**2+p['M_N']*p['g']*p['v_yT1c']*p['x_NG']*sin(q[1])-p['M_N']*p['g']*p['v_yT1c']*p['z_NG']*cos(q[1])-p['M_N']*p['g']*cos(q[1])+p['M_N']*p['v_yT1c']*p['x_NG']*q[2]*qd[1]**2-p['M_N']*p['v_yT1c']*p['z_NG']*(-qd[1]*qd[2]-qd[1]*(p['L_T']*qd[1]+qd[2]))+2*p['M_N']*p['v_yT1c']*q[2]*(p['v_yT1c']*qd[2]+qd[1])*(-p['v_yT1c']*p['x_NG']*qd[2]-p['x_NG']*qd[1])-p['M_N']*qd[1]**2*(-p['v_yT1c']*p['x_NG']*q[2]+p['z_NG'])-p['M_N']*(p['L_T']*p['v_yT1c']-1)*(p['v_yT1c']*qd[2]+qd[1])*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])-p['M_N']*(-qd[1]*qd[2]-qd[1]*(p['L_T']*qd[1]+qd[2]))+p['M_R']*p['g']*p['v_yT1c']*p['x_NR']*sin(q[1])-p['M_R']*p['g']*p['v_yT1c']*p['z_NR']*cos(q[1])-p['M_R']*p['g']*cos(q[1])+p['M_R']*p['v_yT1c']*p['x_NR']*q[2]*qd[1]**2-2*p['M_R']*p['v_yT1c']*p['x_NR']*q[2]*(p['v_yT1c']*qd[2]+qd[1])**2-p['M_R']*p['v_yT1c']*p['z_NR']*(-qd[1]*qd[2]-qd[1]*(p['L_T']*qd[1]+qd[2]))-p['M_R']*p['z_NR']*(p['L_T']*p['v_yT1c']-1)*(p['v_yT1c']*qd[2]+qd[1])**2-p['M_R']*qd[1]**2*(-p['v_yT1c']*p['x_NR']*q[2]+p['z_NR'])-p['M_R']*(-qd[1]*qd[2]-qd[1]*(p['L_T']*qd[1]+qd[2]))
    KK[1,3] = 0
    KK[2,0] = 0
    KK[2,1] = -p['MM_T'][0,6]*p['g']*cos(q[1])-p['M_N']*p['g']*p['v_yT1c']*p['x_NG']*(-p['v_yT1c']*q[2]*cos(q[1])-sin(q[1]))+p['M_N']*p['g']*p['v_yT1c']*p['z_NG']*(p['v_yT1c']*q[2]*sin(q[1])-cos(q[1]))-p['M_N']*p['g']*cos(q[1])-p['M_R']*p['g']*p['v_yT1c']*p['x_NR']*(-p['v_yT1c']*q[2]*cos(q[1])-sin(q[1]))+p['M_R']*p['g']*p['v_yT1c']*p['z_NR']*(p['v_yT1c']*q[2]*sin(q[1])-cos(q[1]))-p['M_R']*p['g']*cos(q[1])
    KK[2,2] = p['KK_T'][6,6]+p['M_N']*p['g']*p['v_yT1c']**2*p['x_NG']*sin(q[1])-p['M_N']*p['g']*p['v_yT1c']**2*p['z_NG']*cos(q[1])+p['M_N']*p['v_yT1c']**2*p['x_NG']*q[2]*qd[1]**2-p['M_N']*p['v_yT1c']**2*p['z_NG']*(-qd[1]*qd[2]-qd[1]*(p['L_T']*qd[1]+qd[2]))-p['M_N']*p['v_yT1c']*(p['v_yT1c']*qd[2]+qd[1])*(p['v_yT1c']*p['z_NG']*qd[2]+p['z_NG']*qd[1])-p['M_N']*qd[1]**2*(-p['v_yT1c']**2*p['x_NG']*q[2]+p['v_yT1c']*p['z_NG'])-p['M_N']*qd[1]**2+p['M_R']*p['g']*p['v_yT1c']**2*p['x_NR']*sin(q[1])-p['M_R']*p['g']*p['v_yT1c']**2*p['z_NR']*cos(q[1])+p['M_R']*p['v_yT1c']**2*p['x_NR']*q[2]*qd[1]**2-p['M_R']*p['v_yT1c']**2*p['z_NR']*(-qd[1]*qd[2]-qd[1]*(p['L_T']*qd[1]+qd[2]))-p['M_R']*p['v_yT1c']*p['z_NR']*(p['v_yT1c']*qd[2]+qd[1])**2-p['M_R']*qd[1]**2*(-p['v_yT1c']**2*p['x_NR']*q[2]+p['v_yT1c']*p['z_NR'])-p['M_R']*qd[1]**2
    KK[2,3] = 0
    KK[3,0] = 0
    KK[3,1] = 0
    KK[3,2] = 0
    KK[3,3] = 0
    return KK

def B_lin(q=None,qd=None,p=None,u=None):
    """ Linear mass matrix 
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'phi_y(t)', 'q_T1(t)', 'psi(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: ['z_T0']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    The columns of B correspond to:   [F_hx(t), F_hz(t)]\\ 
    """
    BB = np.zeros((4,2))
    BB[0,0] = 1
    BB[0,1] = 0
    BB[1,0] = p['z_T0']*cos(q[1])
    BB[1,1] = -p['z_T0']*sin(q[1])
    BB[2,0] = 0
    BB[2,1] = 0
    BB[3,0] = 0
    BB[3,1] = 0
    return BB

