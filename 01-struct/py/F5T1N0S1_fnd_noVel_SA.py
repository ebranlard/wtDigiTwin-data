"""
<YAMSModel object "F5T1N0S1_fnd" with attributes:>
 - coordinates:       [x(t), y(t), phi_x(t), phi_y(t), phi_z(t), q_T1(t), psi(t)]
 - speeds:            [xd(t), yd(t), omega_x_T(t), omega_y_T(t), omega_z_T(t), qd_T1(t), omega_x_R(t)]
 - kdeqsSubs:         [(xd(t), Derivative(x(t), t)), (yd(t), Derivative(y(t), t)), (omega_x_T(t), Derivative(phi_x(t), t)), (omega_y_T(t), Derivative(phi_y(t), t)), (omega_z_T(t), Derivative(phi_z(t), t)), (qd_T1(t), Derivative(q_T1(t), t)), (omega_x_R(t), Derivative(psi(t), t))]
 - var:               [T_a(t)]
 - smallAngleUsed   : [phi_x(t), phi_y(t), phi_z(t)]
 - number of bodies : 4
 - opts             : {'rot_elastic_type': 'SmallRot', 'rot_elastic_smallAngle': False, 'orderMM': 1, 'orderH': 1, 'fnd_loads': False, 'aero_torques': False, 'mergeFndTwr': False, 'yaw': 'zero', 'tilt': 'fixed', 'tiltShaft': True, 'linRot': True, 'Mform': 'symbolic', 'twrDOFDir': ['x', 'y', 'x', 'y'], 'floating': True, 'rot_elastic_subs': True, 'aero_forces': True, 'verbose': False}
 * loads            : [(G_F, - M_F*g*e_E.z), (G_T, - M_T*g*e_E.z), (G_N, - M_N*g*e_E.z), (G_R, - M_R*g*e_E.z), (O_R, T_a(t)*e_R.x)]

"""
import numpy as np
from numpy import cos, sin
def info():
    """ Return information about current model present in this package """
    I=dict()
    I['name']='F5T1N0S1_fnd'
    I['nq']=7
    I['nu']=1
    I['sq']=['x','y','phi_x','phi_y','phi_z','q_T1','psi']
    I['su']=['T_a']
    return I

def M_lin_sa(q=None,p=None,z=None):
    """ Linear mass matrix with small angle approximation
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'y(t)', 'phi_x(t)', 'phi_y(t)', 'phi_z(t)', 'q_T1(t)', 'psi(t)']
    p:  parameters, dictionary with keys: ['JO_R', 'J_xx_F', 'J_xx_N', 'J_yy_F', 'J_yy_N', 'J_zx_N', 'J_zz_F', 'J_zz_N', 'Jxx_R', 'L_T', 'MM_T', 'M_F', 'M_N', 'M_R', 'tilt', 'v_yT1c', 'x_NG', 'x_NR', 'z_FG', 'z_NG', 'z_NR']
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
    MM = np.zeros((7,7))
    MM[0,0] = p['M_F']+p['M_N']+p['M_R']+p['MM_T'][0,0]
    MM[0,1] = 0
    MM[0,2] = 0
    MM[0,3] = p['M_F']*p['z_FG']+p['M_N']*(p['L_T']-p['v_yT1c']*p['x_NG']*q[5]+p['z_NG'])+p['M_R']*(p['L_T']-p['v_yT1c']*p['x_NR']*q[5]+p['z_NR'])-p['MM_T'][1,3]-(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5]))*q[3]
    MM[0,4] = -(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5]))*q[4]
    MM[0,5] = p['MM_T'][0,6]+p['M_N']*(-p['v_yT1c']**2*p['x_NG']*q[5]+p['v_yT1c']*p['z_NG']+1)+p['M_R']*(-p['v_yT1c']**2*p['x_NR']*q[5]+p['v_yT1c']*p['z_NR']+1)-p['v_yT1c']*(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG'])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']))*q[3]
    MM[0,6] = 0
    MM[1,0] = 0
    MM[1,1] = p['M_F']+p['M_N']+p['M_R']+p['MM_T'][0,0]
    MM[1,2] = -p['M_F']*p['z_FG']-p['M_N']*(p['L_T']-p['v_yT1c']*p['x_NG']*q[5]+p['z_NG'])-p['M_R']*(p['L_T']-p['v_yT1c']*p['x_NR']*q[5]+p['z_NR'])+p['MM_T'][1,3]+(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5]))*q[3]
    MM[1,3] = (p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5]))*q[2]
    MM[1,4] = p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5])
    MM[1,5] = p['v_yT1c']*(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG'])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']))*q[2]+(p['MM_T'][0,6]+p['M_N']*(-p['v_yT1c']**2*p['x_NG']*q[5]+p['v_yT1c']*p['z_NG']+1)+p['M_R']*(-p['v_yT1c']**2*p['x_NR']*q[5]+p['v_yT1c']*p['z_NR']+1))*q[4]
    MM[1,6] = 0
    MM[2,0] = 0
    MM[2,1] = -p['M_F']*p['z_FG']-p['M_N']*(p['L_T']-p['v_yT1c']*p['x_NG']*q[5]+p['z_NG'])-p['M_R']*(p['L_T']-p['v_yT1c']*p['x_NR']*q[5]+p['z_NR'])+p['MM_T'][1,3]+(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5]))*q[3]
    MM[2,2] = -p['JO_R']*p['v_yT1c']**2*q[3]*q[5]**2*sin(2*p['tilt'])+p['JO_R']*p['v_yT1c']**2*q[5]**2*cos(2*p['tilt'])/2+p['JO_R']*p['v_yT1c']**2*q[5]**2/2+2*p['JO_R']*p['v_yT1c']*q[3]*q[5]*cos(2*p['tilt'])+p['JO_R']*p['v_yT1c']*q[5]*sin(2*p['tilt'])+p['JO_R']*q[3]*sin(2*p['tilt'])-p['JO_R']*cos(2*p['tilt'])/2+p['JO_R']/2+p['MM_T'][3,3]+p['J_xx_F']-2*p['J_xx_N']*p['v_yT1c']*q[3]*q[5]+p['J_xx_N']-2*p['J_zx_N']*p['v_yT1c']**2*q[3]*q[5]**2+2*p['J_zx_N']*p['v_yT1c']*q[5]+2*p['J_zx_N']*q[3]+p['J_zz_N']*p['v_yT1c']**2*q[5]**2+2*p['J_zz_N']*p['v_yT1c']*q[3]*q[5]+p['Jxx_R']*p['v_yT1c']**2*q[3]*q[5]**2*sin(2*p['tilt'])-p['Jxx_R']*p['v_yT1c']**2*q[5]**2*cos(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']**2*q[5]**2/2-2*p['Jxx_R']*p['v_yT1c']*q[3]*q[5]*cos(2*p['tilt'])-p['Jxx_R']*p['v_yT1c']*q[5]*sin(2*p['tilt'])-p['Jxx_R']*q[3]*sin(2*p['tilt'])+p['Jxx_R']*cos(2*p['tilt'])/2+p['Jxx_R']/2+p['L_T']**2*p['M_N']+p['L_T']**2*p['M_R']-2*p['L_T']*p['M_N']*p['v_yT1c']*p['x_NG']*q[5]-2*p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']*q[3]*q[5]-2*p['L_T']*p['M_N']*p['x_NG']*q[3]+2*p['L_T']*p['M_N']*p['z_NG']-2*p['L_T']*p['M_N']*q[3]*q[5]-2*p['L_T']*p['M_R']*p['v_yT1c']*p['x_NR']*q[5]-2*p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']*q[3]*q[5]-2*p['L_T']*p['M_R']*p['x_NR']*q[3]+2*p['L_T']*p['M_R']*p['z_NR']-2*p['L_T']*p['M_R']*q[3]*q[5]+p['M_F']*p['z_FG']**2+p['M_N']*p['v_yT1c']**2*p['x_NG']**2*q[5]**2+2*p['M_N']*p['v_yT1c']**2*p['x_NG']*p['z_NG']*q[3]*q[5]**2+2*p['M_N']*p['v_yT1c']*p['x_NG']**2*q[3]*q[5]-2*p['M_N']*p['v_yT1c']*p['x_NG']*p['z_NG']*q[5]+2*p['M_N']*p['v_yT1c']*p['x_NG']*q[3]*q[5]**2-2*p['M_N']*p['v_yT1c']*p['z_NG']**2*q[3]*q[5]-2*p['M_N']*p['x_NG']*p['z_NG']*q[3]+p['M_N']*p['z_NG']**2-2*p['M_N']*p['z_NG']*q[3]*q[5]+p['M_R']*p['v_yT1c']**2*p['x_NR']**2*q[5]**2+2*p['M_R']*p['v_yT1c']**2*p['x_NR']*p['z_NR']*q[3]*q[5]**2+2*p['M_R']*p['v_yT1c']*p['x_NR']**2*q[3]*q[5]-2*p['M_R']*p['v_yT1c']*p['x_NR']*p['z_NR']*q[5]+2*p['M_R']*p['v_yT1c']*p['x_NR']*q[3]*q[5]**2-2*p['M_R']*p['v_yT1c']*p['z_NR']**2*q[3]*q[5]-2*p['M_R']*p['x_NR']*p['z_NR']*q[3]+p['M_R']*p['z_NR']**2-2*p['M_R']*p['z_NR']*q[3]*q[5]
    MM[2,3] = (p['JO_R']*p['v_yT1c']**2*q[5]**2*cos(2*p['tilt'])/2+p['JO_R']*p['v_yT1c']**2*q[5]**2/2+p['JO_R']*p['v_yT1c']*q[5]*sin(2*p['tilt'])-p['JO_R']*cos(2*p['tilt'])/2-p['JO_R']/2+p['MM_T'][3,3]-p['MM_T'][4,4]+p['J_xx_F']+p['J_xx_N']-p['J_yy_F']-p['J_yy_N']+2*p['J_zx_N']*p['v_yT1c']*q[5]+p['J_zz_N']*p['v_yT1c']**2*q[5]**2-p['Jxx_R']*p['v_yT1c']**2*q[5]**2*cos(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']**2*q[5]**2/2-p['Jxx_R']*p['v_yT1c']*q[5]*sin(2*p['tilt'])+p['Jxx_R']*cos(2*p['tilt'])/2+p['Jxx_R']/2+p['M_N']*p['v_yT1c']**2*p['x_NG']**2*q[5]**2-2*p['M_N']*p['v_yT1c']*p['x_NG']*p['z_NG']*q[5]-2*p['M_N']*p['v_yT1c']*p['z_NG']*q[5]**2-p['M_N']*p['x_NG']**2-2*p['M_N']*p['x_NG']*q[5]-p['M_N']*q[5]**2+p['M_R']*p['v_yT1c']**2*p['x_NR']**2*q[5]**2-2*p['M_R']*p['v_yT1c']*p['x_NR']*p['z_NR']*q[5]-2*p['M_R']*p['v_yT1c']*p['z_NR']*q[5]**2-p['M_R']*p['x_NR']**2-2*p['M_R']*p['x_NR']*q[5]-p['M_R']*q[5]**2)*q[4]
    MM[2,4] = -p['JO_R']*p['v_yT1c']**2*q[3]*q[5]**2*cos(2*p['tilt'])/2+p['JO_R']*p['v_yT1c']**2*q[3]*q[5]**2/2-p['JO_R']*p['v_yT1c']**2*q[5]**2*sin(2*p['tilt'])/2-p['JO_R']*p['v_yT1c']*q[3]*q[5]*sin(2*p['tilt'])+p['JO_R']*p['v_yT1c']*q[5]*cos(2*p['tilt'])+p['JO_R']*q[3]*cos(2*p['tilt'])/2+p['JO_R']*q[3]/2+p['JO_R']*sin(2*p['tilt'])/2+p['MM_T'][5,5]*q[3]+p['J_xx_N']*p['v_yT1c']**2*q[3]*q[5]**2-p['J_xx_N']*p['v_yT1c']*q[5]-p['J_zx_N']*p['v_yT1c']**2*q[5]**2-2*p['J_zx_N']*p['v_yT1c']*q[3]*q[5]+p['J_zx_N']+p['J_zz_F']*q[3]+p['J_zz_N']*p['v_yT1c']*q[5]+p['J_zz_N']*q[3]+p['Jxx_R']*p['v_yT1c']**2*q[3]*q[5]**2*cos(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']**2*q[3]*q[5]**2/2+p['Jxx_R']*p['v_yT1c']**2*q[5]**2*sin(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']*q[3]*q[5]*sin(2*p['tilt'])-p['Jxx_R']*p['v_yT1c']*q[5]*cos(2*p['tilt'])-p['Jxx_R']*q[3]*cos(2*p['tilt'])/2+p['Jxx_R']*q[3]/2-p['Jxx_R']*sin(2*p['tilt'])/2-p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']*q[5]-p['L_T']*p['M_N']*p['x_NG']-p['L_T']*p['M_N']*q[5]-p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']*q[5]-p['L_T']*p['M_R']*p['x_NR']-p['L_T']*p['M_R']*q[5]+p['M_N']*p['v_yT1c']**2*p['x_NG']*p['z_NG']*q[5]**2+p['M_N']*p['v_yT1c']**2*p['z_NG']**2*q[3]*q[5]**2+p['M_N']*p['v_yT1c']*p['x_NG']**2*q[5]+2*p['M_N']*p['v_yT1c']*p['x_NG']*p['z_NG']*q[3]*q[5]+p['M_N']*p['v_yT1c']*p['x_NG']*q[5]**2-p['M_N']*p['v_yT1c']*p['z_NG']**2*q[5]+2*p['M_N']*p['v_yT1c']*p['z_NG']*q[3]*q[5]**2+p['M_N']*p['x_NG']**2*q[3]-p['M_N']*p['x_NG']*p['z_NG']+2*p['M_N']*p['x_NG']*q[3]*q[5]-p['M_N']*p['z_NG']*q[5]+p['M_N']*q[3]*q[5]**2+p['M_R']*p['v_yT1c']**2*p['x_NR']*p['z_NR']*q[5]**2+p['M_R']*p['v_yT1c']**2*p['z_NR']**2*q[3]*q[5]**2+p['M_R']*p['v_yT1c']*p['x_NR']**2*q[5]+2*p['M_R']*p['v_yT1c']*p['x_NR']*p['z_NR']*q[3]*q[5]+p['M_R']*p['v_yT1c']*p['x_NR']*q[5]**2-p['M_R']*p['v_yT1c']*p['z_NR']**2*q[5]+2*p['M_R']*p['v_yT1c']*p['z_NR']*q[3]*q[5]**2+p['M_R']*p['x_NR']**2*q[3]-p['M_R']*p['x_NR']*p['z_NR']+2*p['M_R']*p['x_NR']*q[3]*q[5]-p['M_R']*p['z_NR']*q[5]+p['M_R']*q[3]*q[5]**2
    MM[2,5] = -(p['MM_T'][6,4]+p['JO_R']*p['v_yT1c']+p['J_yy_N']*p['v_yT1c']-p['L_T']*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[5]+p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']+p['L_T']*p['M_N']-p['L_T']*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[5]+p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']+p['L_T']*p['M_R']+p['M_N']*p['v_yT1c']**2*p['z_NG']*q[5]**2+p['M_N']*p['v_yT1c']*p['x_NG']**2+p['M_N']*p['v_yT1c']*p['z_NG']**2+p['M_N']*p['z_NG']+p['M_R']*p['v_yT1c']**2*p['z_NR']*q[5]**2+p['M_R']*p['v_yT1c']*p['x_NR']**2+p['M_R']*p['v_yT1c']*p['z_NR']**2+p['M_R']*p['z_NR'])*q[4]
    MM[2,6] = p['Jxx_R']*(-p['v_yT1c']*q[5]*sin(p['tilt'])-(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*q[3]+cos(p['tilt']))
    MM[3,0] = p['M_F']*p['z_FG']+p['M_N']*(p['L_T']-p['v_yT1c']*p['x_NG']*q[5]+p['z_NG'])+p['M_R']*(p['L_T']-p['v_yT1c']*p['x_NR']*q[5]+p['z_NR'])-p['MM_T'][1,3]-(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5]))*q[3]
    MM[3,1] = (p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5]))*q[2]
    MM[3,2] = (p['JO_R']*p['v_yT1c']**2*q[5]**2*cos(2*p['tilt'])/2+p['JO_R']*p['v_yT1c']**2*q[5]**2/2+p['JO_R']*p['v_yT1c']*q[5]*sin(2*p['tilt'])-p['JO_R']*cos(2*p['tilt'])/2-p['JO_R']/2+p['MM_T'][3,3]-p['MM_T'][4,4]+p['J_xx_F']+p['J_xx_N']-p['J_yy_F']-p['J_yy_N']+2*p['J_zx_N']*p['v_yT1c']*q[5]+p['J_zz_N']*p['v_yT1c']**2*q[5]**2-p['Jxx_R']*p['v_yT1c']**2*q[5]**2*cos(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']**2*q[5]**2/2-p['Jxx_R']*p['v_yT1c']*q[5]*sin(2*p['tilt'])+p['Jxx_R']*cos(2*p['tilt'])/2+p['Jxx_R']/2+p['M_N']*p['v_yT1c']**2*p['x_NG']**2*q[5]**2-2*p['M_N']*p['v_yT1c']*p['x_NG']*p['z_NG']*q[5]-2*p['M_N']*p['v_yT1c']*p['z_NG']*q[5]**2-p['M_N']*p['x_NG']**2-2*p['M_N']*p['x_NG']*q[5]-p['M_N']*q[5]**2+p['M_R']*p['v_yT1c']**2*p['x_NR']**2*q[5]**2-2*p['M_R']*p['v_yT1c']*p['x_NR']*p['z_NR']*q[5]-2*p['M_R']*p['v_yT1c']*p['z_NR']*q[5]**2-p['M_R']*p['x_NR']**2-2*p['M_R']*p['x_NR']*q[5]-p['M_R']*q[5]**2)*q[4]
    MM[3,3] = p['JO_R']+p['MM_T'][4,4]+p['J_yy_F']+p['J_yy_N']+p['L_T']**2*p['M_N']+p['L_T']**2*p['M_R']-2*p['L_T']*p['M_N']*p['v_yT1c']*p['x_NG']*q[5]+2*p['L_T']*p['M_N']*p['z_NG']-2*p['L_T']*p['M_R']*p['v_yT1c']*p['x_NR']*q[5]+2*p['L_T']*p['M_R']*p['z_NR']+p['M_F']*p['z_FG']**2+2*p['M_N']*p['v_yT1c']*p['z_NG']*q[5]**2+p['M_N']*p['x_NG']**2+2*p['M_N']*p['x_NG']*q[5]+p['M_N']*p['z_NG']**2+p['M_N']*q[5]**2+2*p['M_R']*p['v_yT1c']*p['z_NR']*q[5]**2+p['M_R']*p['x_NR']**2+2*p['M_R']*p['x_NR']*q[5]+p['M_R']*p['z_NR']**2+p['M_R']*q[5]**2
    MM[3,4] = -(p['JO_R']*p['v_yT1c']**2*q[5]**2*sin(2*p['tilt'])/2-p['JO_R']*p['v_yT1c']*q[5]*cos(2*p['tilt'])-p['JO_R']*sin(2*p['tilt'])/2+p['J_xx_N']*p['v_yT1c']*q[5]+p['J_zx_N']*p['v_yT1c']**2*q[5]**2-p['J_zx_N']-p['J_zz_N']*p['v_yT1c']*q[5]-p['Jxx_R']*p['v_yT1c']**2*q[5]**2*sin(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']*q[5]*cos(2*p['tilt'])+p['Jxx_R']*sin(2*p['tilt'])/2+p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']*q[5]+p['L_T']*p['M_N']*p['x_NG']+p['L_T']*p['M_N']*q[5]+p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']*q[5]+p['L_T']*p['M_R']*p['x_NR']+p['L_T']*p['M_R']*q[5]-p['M_N']*p['v_yT1c']**2*p['x_NG']*p['z_NG']*q[5]**2-p['M_N']*p['v_yT1c']*p['x_NG']**2*q[5]-p['M_N']*p['v_yT1c']*p['x_NG']*q[5]**2+p['M_N']*p['v_yT1c']*p['z_NG']**2*q[5]+p['M_N']*p['x_NG']*p['z_NG']+p['M_N']*p['z_NG']*q[5]-p['M_R']*p['v_yT1c']**2*p['x_NR']*p['z_NR']*q[5]**2-p['M_R']*p['v_yT1c']*p['x_NR']**2*q[5]-p['M_R']*p['v_yT1c']*p['x_NR']*q[5]**2+p['M_R']*p['v_yT1c']*p['z_NR']**2*q[5]+p['M_R']*p['x_NR']*p['z_NR']+p['M_R']*p['z_NR']*q[5])*q[4]
    MM[3,5] = p['MM_T'][6,4]+p['JO_R']*p['v_yT1c']+p['J_yy_N']*p['v_yT1c']-p['L_T']*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[5]+p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']+p['L_T']*p['M_N']-p['L_T']*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[5]+p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']+p['L_T']*p['M_R']+p['M_N']*p['v_yT1c']**2*p['z_NG']*q[5]**2+p['M_N']*p['v_yT1c']*p['x_NG']**2+p['M_N']*p['v_yT1c']*p['z_NG']**2+p['M_N']*p['z_NG']+p['M_R']*p['v_yT1c']**2*p['z_NR']*q[5]**2+p['M_R']*p['v_yT1c']*p['x_NR']**2+p['M_R']*p['v_yT1c']*p['z_NR']**2+p['M_R']*p['z_NR']
    MM[3,6] = p['Jxx_R']*(-p['v_yT1c']*q[5]*sin(p['tilt'])+cos(p['tilt']))*q[4]
    MM[4,0] = -(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5]))*q[4]
    MM[4,1] = p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG']+q[5])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']+q[5])
    MM[4,2] = -p['JO_R']*p['v_yT1c']**2*q[3]*q[5]**2*cos(2*p['tilt'])/2+p['JO_R']*p['v_yT1c']**2*q[3]*q[5]**2/2-p['JO_R']*p['v_yT1c']**2*q[5]**2*sin(2*p['tilt'])/2-p['JO_R']*p['v_yT1c']*q[3]*q[5]*sin(2*p['tilt'])+p['JO_R']*p['v_yT1c']*q[5]*cos(2*p['tilt'])+p['JO_R']*q[3]*cos(2*p['tilt'])/2+p['JO_R']*q[3]/2+p['JO_R']*sin(2*p['tilt'])/2+p['MM_T'][5,5]*q[3]+p['J_xx_N']*p['v_yT1c']**2*q[3]*q[5]**2-p['J_xx_N']*p['v_yT1c']*q[5]-p['J_zx_N']*p['v_yT1c']**2*q[5]**2-2*p['J_zx_N']*p['v_yT1c']*q[3]*q[5]+p['J_zx_N']+p['J_zz_F']*q[3]+p['J_zz_N']*p['v_yT1c']*q[5]+p['J_zz_N']*q[3]+p['Jxx_R']*p['v_yT1c']**2*q[3]*q[5]**2*cos(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']**2*q[3]*q[5]**2/2+p['Jxx_R']*p['v_yT1c']**2*q[5]**2*sin(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']*q[3]*q[5]*sin(2*p['tilt'])-p['Jxx_R']*p['v_yT1c']*q[5]*cos(2*p['tilt'])-p['Jxx_R']*q[3]*cos(2*p['tilt'])/2+p['Jxx_R']*q[3]/2-p['Jxx_R']*sin(2*p['tilt'])/2-p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']*q[5]-p['L_T']*p['M_N']*p['x_NG']-p['L_T']*p['M_N']*q[5]-p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']*q[5]-p['L_T']*p['M_R']*p['x_NR']-p['L_T']*p['M_R']*q[5]+p['M_N']*p['v_yT1c']**2*p['x_NG']*p['z_NG']*q[5]**2+p['M_N']*p['v_yT1c']**2*p['z_NG']**2*q[3]*q[5]**2+p['M_N']*p['v_yT1c']*p['x_NG']**2*q[5]+2*p['M_N']*p['v_yT1c']*p['x_NG']*p['z_NG']*q[3]*q[5]+p['M_N']*p['v_yT1c']*p['x_NG']*q[5]**2-p['M_N']*p['v_yT1c']*p['z_NG']**2*q[5]+2*p['M_N']*p['v_yT1c']*p['z_NG']*q[3]*q[5]**2+p['M_N']*p['x_NG']**2*q[3]-p['M_N']*p['x_NG']*p['z_NG']+2*p['M_N']*p['x_NG']*q[3]*q[5]-p['M_N']*p['z_NG']*q[5]+p['M_N']*q[3]*q[5]**2+p['M_R']*p['v_yT1c']**2*p['x_NR']*p['z_NR']*q[5]**2+p['M_R']*p['v_yT1c']**2*p['z_NR']**2*q[3]*q[5]**2+p['M_R']*p['v_yT1c']*p['x_NR']**2*q[5]+2*p['M_R']*p['v_yT1c']*p['x_NR']*p['z_NR']*q[3]*q[5]+p['M_R']*p['v_yT1c']*p['x_NR']*q[5]**2-p['M_R']*p['v_yT1c']*p['z_NR']**2*q[5]+2*p['M_R']*p['v_yT1c']*p['z_NR']*q[3]*q[5]**2+p['M_R']*p['x_NR']**2*q[3]-p['M_R']*p['x_NR']*p['z_NR']+2*p['M_R']*p['x_NR']*q[3]*q[5]-p['M_R']*p['z_NR']*q[5]+p['M_R']*q[3]*q[5]**2
    MM[4,3] = -(p['JO_R']*p['v_yT1c']**2*q[5]**2*sin(2*p['tilt'])/2-p['JO_R']*p['v_yT1c']*q[5]*cos(2*p['tilt'])-p['JO_R']*sin(2*p['tilt'])/2+p['J_xx_N']*p['v_yT1c']*q[5]+p['J_zx_N']*p['v_yT1c']**2*q[5]**2-p['J_zx_N']-p['J_zz_N']*p['v_yT1c']*q[5]-p['Jxx_R']*p['v_yT1c']**2*q[5]**2*sin(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']*q[5]*cos(2*p['tilt'])+p['Jxx_R']*sin(2*p['tilt'])/2+p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']*q[5]+p['L_T']*p['M_N']*p['x_NG']+p['L_T']*p['M_N']*q[5]+p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']*q[5]+p['L_T']*p['M_R']*p['x_NR']+p['L_T']*p['M_R']*q[5]-p['M_N']*p['v_yT1c']**2*p['x_NG']*p['z_NG']*q[5]**2-p['M_N']*p['v_yT1c']*p['x_NG']**2*q[5]-p['M_N']*p['v_yT1c']*p['x_NG']*q[5]**2+p['M_N']*p['v_yT1c']*p['z_NG']**2*q[5]+p['M_N']*p['x_NG']*p['z_NG']+p['M_N']*p['z_NG']*q[5]-p['M_R']*p['v_yT1c']**2*p['x_NR']*p['z_NR']*q[5]**2-p['M_R']*p['v_yT1c']*p['x_NR']**2*q[5]-p['M_R']*p['v_yT1c']*p['x_NR']*q[5]**2+p['M_R']*p['v_yT1c']*p['z_NR']**2*q[5]+p['M_R']*p['x_NR']*p['z_NR']+p['M_R']*p['z_NR']*q[5])*q[4]
    MM[4,4] = -p['JO_R']*p['v_yT1c']**2*q[5]**2*cos(2*p['tilt'])/2+p['JO_R']*p['v_yT1c']**2*q[5]**2/2-p['JO_R']*p['v_yT1c']*q[5]*sin(2*p['tilt'])+p['JO_R']*cos(2*p['tilt'])/2+p['JO_R']/2+p['MM_T'][5,5]+p['J_xx_N']*p['v_yT1c']**2*q[5]**2-2*p['J_zx_N']*p['v_yT1c']*q[5]+p['J_zz_F']+p['J_zz_N']+p['Jxx_R']*p['v_yT1c']**2*q[5]**2*cos(2*p['tilt'])/2+p['Jxx_R']*p['v_yT1c']**2*q[5]**2/2+p['Jxx_R']*p['v_yT1c']*q[5]*sin(2*p['tilt'])-p['Jxx_R']*cos(2*p['tilt'])/2+p['Jxx_R']/2+p['M_N']*p['v_yT1c']**2*p['z_NG']**2*q[5]**2+2*p['M_N']*p['v_yT1c']*p['x_NG']*p['z_NG']*q[5]+2*p['M_N']*p['v_yT1c']*p['z_NG']*q[5]**2+p['M_N']*p['x_NG']**2+2*p['M_N']*p['x_NG']*q[5]+p['M_N']*q[5]**2+p['M_R']*p['v_yT1c']**2*p['z_NR']**2*q[5]**2+2*p['M_R']*p['v_yT1c']*p['x_NR']*p['z_NR']*q[5]+2*p['M_R']*p['v_yT1c']*p['z_NR']*q[5]**2+p['M_R']*p['x_NR']**2+2*p['M_R']*p['x_NR']*q[5]+p['M_R']*q[5]**2
    MM[4,5] = 0
    MM[4,6] = -p['Jxx_R']*(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))
    MM[5,0] = p['MM_T'][0,6]+p['M_N']*(-p['v_yT1c']**2*p['x_NG']*q[5]+p['v_yT1c']*p['z_NG']+1)+p['M_R']*(-p['v_yT1c']**2*p['x_NR']*q[5]+p['v_yT1c']*p['z_NR']+1)-p['v_yT1c']*(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG'])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']))*q[3]
    MM[5,1] = p['v_yT1c']*(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG'])+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR']))*q[2]+(p['MM_T'][0,6]+p['M_N']*(-p['v_yT1c']**2*p['x_NG']*q[5]+p['v_yT1c']*p['z_NG']+1)+p['M_R']*(-p['v_yT1c']**2*p['x_NR']*q[5]+p['v_yT1c']*p['z_NR']+1))*q[4]
    MM[5,2] = -(p['MM_T'][6,4]+p['JO_R']*p['v_yT1c']+p['J_yy_N']*p['v_yT1c']-p['L_T']*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[5]+p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']+p['L_T']*p['M_N']-p['L_T']*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[5]+p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']+p['L_T']*p['M_R']+p['M_N']*p['v_yT1c']**2*p['z_NG']*q[5]**2+p['M_N']*p['v_yT1c']*p['x_NG']**2+p['M_N']*p['v_yT1c']*p['z_NG']**2+p['M_N']*p['z_NG']+p['M_R']*p['v_yT1c']**2*p['z_NR']*q[5]**2+p['M_R']*p['v_yT1c']*p['x_NR']**2+p['M_R']*p['v_yT1c']*p['z_NR']**2+p['M_R']*p['z_NR'])*q[4]
    MM[5,3] = p['MM_T'][6,4]+p['JO_R']*p['v_yT1c']+p['J_yy_N']*p['v_yT1c']-p['L_T']*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[5]+p['L_T']*p['M_N']*p['v_yT1c']*p['z_NG']+p['L_T']*p['M_N']-p['L_T']*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[5]+p['L_T']*p['M_R']*p['v_yT1c']*p['z_NR']+p['L_T']*p['M_R']+p['M_N']*p['v_yT1c']**2*p['z_NG']*q[5]**2+p['M_N']*p['v_yT1c']*p['x_NG']**2+p['M_N']*p['v_yT1c']*p['z_NG']**2+p['M_N']*p['z_NG']+p['M_R']*p['v_yT1c']**2*p['z_NR']*q[5]**2+p['M_R']*p['v_yT1c']*p['x_NR']**2+p['M_R']*p['v_yT1c']*p['z_NR']**2+p['M_R']*p['z_NR']
    MM[5,4] = 0
    MM[5,5] = p['JO_R']*p['v_yT1c']**2+p['J_yy_N']*p['v_yT1c']**2+p['M_N']*p['v_yT1c']**2*p['x_NG']**2-2*p['M_N']*p['v_yT1c']**2*p['x_NG']*q[5]+p['M_N']*p['v_yT1c']**2*p['z_NG']**2+2*p['M_N']*p['v_yT1c']*p['z_NG']+p['M_N']+p['M_R']*p['v_yT1c']**2*p['x_NR']**2-2*p['M_R']*p['v_yT1c']**2*p['x_NR']*q[5]+p['M_R']*p['v_yT1c']**2*p['z_NR']**2+2*p['M_R']*p['v_yT1c']*p['z_NR']+p['M_R']+p['MM_T'][6,6]
    MM[5,6] = 0
    MM[6,0] = 0
    MM[6,1] = 0
    MM[6,2] = p['Jxx_R']*(-p['v_yT1c']*q[5]*sin(p['tilt'])-(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*q[3]+cos(p['tilt']))
    MM[6,3] = p['Jxx_R']*(-p['v_yT1c']*q[5]*sin(p['tilt'])+cos(p['tilt']))*q[4]
    MM[6,4] = -p['Jxx_R']*(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))
    MM[6,5] = 0
    MM[6,6] = p['Jxx_R']
    return MM

def C_lin_sa(q=None,qd=None,p=None,u=None,z=None):
    """ Linear damping matrix with small angle approximation
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'y(t)', 'phi_x(t)', 'phi_y(t)', 'phi_z(t)', 'q_T1(t)', 'psi(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: ['DD_T', 'Jxx_R', 'tilt', 'v_yT1c']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    CC = np.zeros((7,7))
    CC[0,0] = 0
    CC[0,1] = 0
    CC[0,2] = 0
    CC[0,3] = 0
    CC[0,4] = 0
    CC[0,5] = 0
    CC[0,6] = 0
    CC[1,0] = 0
    CC[1,1] = 0
    CC[1,2] = 0
    CC[1,3] = 0
    CC[1,4] = 0
    CC[1,5] = 0
    CC[1,6] = 0
    CC[2,0] = 0
    CC[2,1] = 0
    CC[2,2] = 0
    CC[2,3] = p['Jxx_R']*(p['v_yT1c']*q[3]*q[5]*sin(p['tilt'])-p['v_yT1c']*q[5]*cos(p['tilt'])-q[3]*cos(p['tilt'])-sin(p['tilt']))*qd[6]
    CC[2,4] = p['Jxx_R']*(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))*q[4]*qd[6]
    CC[2,5] = p['Jxx_R']*p['v_yT1c']*(p['v_yT1c']*q[3]*q[5]*sin(p['tilt'])-p['v_yT1c']*q[5]*cos(p['tilt'])-q[3]*cos(p['tilt'])-sin(p['tilt']))*qd[6]
    CC[2,6] = 0
    CC[3,0] = 0
    CC[3,1] = 0
    CC[3,2] = p['Jxx_R']*(-p['v_yT1c']*q[3]*q[5]*sin(p['tilt'])+p['v_yT1c']*q[5]*cos(p['tilt'])+q[3]*cos(p['tilt'])+sin(p['tilt']))*qd[6]
    CC[3,3] = 0
    CC[3,4] = -p['Jxx_R']*(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))*qd[6]
    CC[3,5] = -p['Jxx_R']*p['v_yT1c']*(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*q[4]*qd[6]
    CC[3,6] = 0
    CC[4,0] = 0
    CC[4,1] = 0
    CC[4,2] = -p['Jxx_R']*(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))*q[4]*qd[6]
    CC[4,3] = p['Jxx_R']*(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))*qd[6]
    CC[4,4] = 0
    CC[4,5] = p['Jxx_R']*p['v_yT1c']*(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))*qd[6]
    CC[4,6] = 0
    CC[5,0] = 0
    CC[5,1] = 0
    CC[5,2] = p['Jxx_R']*p['v_yT1c']*(-p['v_yT1c']*q[3]*q[5]*sin(p['tilt'])+p['v_yT1c']*q[5]*cos(p['tilt'])+q[3]*cos(p['tilt'])+sin(p['tilt']))*qd[6]
    CC[5,3] = p['Jxx_R']*p['v_yT1c']*(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*q[4]*qd[6]
    CC[5,4] = -p['Jxx_R']*p['v_yT1c']*(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))*qd[6]
    CC[5,5] = p['DD_T'][6,6]
    CC[5,6] = 0
    CC[6,0] = 0
    CC[6,1] = 0
    CC[6,2] = 0
    CC[6,3] = 0
    CC[6,4] = 0
    CC[6,5] = 0
    CC[6,6] = 0
    return CC

def K_lin_sa(q=None,qd=None,p=None,u=None,z=None):
    """ Linear stiffness matrix with small angle approximation
    q:  degrees of freedom, array-like: ['x(t)', 'y(t)', 'phi_x(t)', 'phi_y(t)', 'phi_z(t)', 'q_T1(t)', 'psi(t)']
    qd: dof velocities, array-like
    p:  parameters, dictionary with keys: ['KK_T', 'L_T', 'MM_T', 'M_F', 'M_N', 'M_R', 'g', 'tilt', 'v_yT1c', 'x_NG', 'x_NR', 'z_FG', 'z_NG', 'z_NR']
    u:  inputs at operating point, dictionary with keys: ['T_a']
           where each values is a constant!
    """
    if z is not None:
        q  = z[0:int(len(z)/2)] 
        qd = z[int(len(z)/2): ] 
    KK = np.zeros((7,7))
    KK[0,0] = 0
    KK[0,1] = 0
    KK[0,2] = 0
    KK[0,3] = (p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*u['T_a']
    KK[0,4] = 0
    KK[0,5] = p['v_yT1c']*(q[3]*cos(p['tilt'])+sin(p['tilt']))*u['T_a']
    KK[0,6] = 0
    KK[1,0] = 0
    KK[1,1] = 0
    KK[1,2] = -(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*u['T_a']
    KK[1,3] = 0
    KK[1,4] = (p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))*u['T_a']
    KK[1,5] = p['v_yT1c']*(-q[2]*cos(p['tilt'])+q[4]*sin(p['tilt']))*u['T_a']
    KK[1,6] = 0
    KK[2,0] = 0
    KK[2,1] = 0
    KK[2,2] = p['g']*(-p['L_T']*p['M_N']-p['L_T']*p['M_R']-p['M_F']*p['z_FG']+p['M_N']*(p['v_yT1c']*p['x_NG']*q[5]-p['z_NG'])+p['M_R']*(p['v_yT1c']*p['x_NR']*q[5]-p['z_NR'])+p['MM_T'][1,3])
    KK[2,3] = 0
    KK[2,4] = -p['L_T']*(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))*u['T_a']+p['M_N']*p['g']*p['v_yT1c']*p['z_NG']*q[5]+p['M_N']*p['g']*p['x_NG']+p['M_N']*p['g']*q[5]+p['M_R']*p['g']*p['v_yT1c']*p['z_NR']*q[5]+p['M_R']*p['g']*p['x_NR']+p['M_R']*p['g']*q[5]+p['x_NR']*u['T_a']*sin(p['tilt'])+p['z_NR']*u['T_a']*cos(p['tilt'])+(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*u['T_a']*q[5]
    KK[2,5] = p['g']*p['v_yT1c']*(p['M_N']*p['x_NG']+p['M_R']*p['x_NR'])*q[2]+(-p['L_T']*p['v_yT1c']*u['T_a']*sin(p['tilt'])+p['M_N']*p['g']*p['v_yT1c']*p['z_NG']+p['M_N']*p['g']+p['M_R']*p['g']*p['v_yT1c']*p['z_NR']+p['M_R']*p['g']+2*p['v_yT1c']*u['T_a']*q[5]*cos(p['tilt'])+u['T_a']*sin(p['tilt']))*q[4]
    KK[2,6] = 0
    KK[3,0] = 0
    KK[3,1] = 0
    KK[3,2] = 0
    KK[3,3] = p['g']*(-p['L_T']*p['M_N']-p['L_T']*p['M_R']-p['M_F']*p['z_FG']+p['M_N']*p['v_yT1c']*p['x_NG']*q[5]-p['M_N']*p['z_NG']+p['M_R']*p['v_yT1c']*p['x_NR']*q[5]-p['M_R']*p['z_NR']+p['MM_T'][1,3])
    KK[3,4] = 0
    KK[3,5] = p['L_T']*p['v_yT1c']*u['T_a']*sin(p['tilt'])+p['M_N']*p['g']*p['v_yT1c']*p['x_NG']*q[3]-p['M_N']*p['g']*p['v_yT1c']*p['z_NG']-p['M_N']*p['g']+p['M_R']*p['g']*p['v_yT1c']*p['x_NR']*q[3]-p['M_R']*p['g']*p['v_yT1c']*p['z_NR']-p['M_R']*p['g']-2*p['v_yT1c']*u['T_a']*q[5]*cos(p['tilt'])-u['T_a']*sin(p['tilt'])
    KK[3,6] = 0
    KK[4,0] = 0
    KK[4,1] = 0
    KK[4,2] = p['g']*(p['M_N']*(p['v_yT1c']*p['z_NG']*q[5]+p['x_NG'])+p['M_N']*q[5]+p['M_R']*(p['v_yT1c']*p['z_NR']*q[5]+p['x_NR'])+p['M_R']*q[5])
    KK[4,3] = 0
    KK[4,4] = 0
    KK[4,5] = p['g']*(p['M_N']*p['v_yT1c']*p['z_NG']+p['M_N']+p['M_R']*p['v_yT1c']*p['z_NR']+p['M_R'])*q[2]
    KK[4,6] = 0
    KK[5,0] = 0
    KK[5,1] = 0
    KK[5,2] = 0
    KK[5,3] = p['g']*(-p['MM_T'][0,6]+p['M_N']*p['v_yT1c']**2*p['x_NG']*q[5]-p['M_N']*p['v_yT1c']*p['z_NG']-p['M_N']+p['M_R']*p['v_yT1c']**2*p['x_NR']*q[5]-p['M_R']*p['v_yT1c']*p['z_NR']-p['M_R'])
    KK[5,4] = 0
    KK[5,5] = p['KK_T'][6,6]-p['M_N']*p['g']*p['v_yT1c']**2*p['z_NG']-p['M_R']*p['g']*p['v_yT1c']**2*p['z_NR']+p['g']*p['v_yT1c']**2*(p['M_N']*p['x_NG']+p['M_R']*p['x_NR'])*q[3]+p['v_yT1c']*u['T_a']*sin(p['tilt'])
    KK[5,6] = 0
    KK[6,0] = 0
    KK[6,1] = 0
    KK[6,2] = 0
    KK[6,3] = 0
    KK[6,4] = 0
    KK[6,5] = 0
    KK[6,6] = 0
    return KK

def B_lin_sa(q=None,qd=None,p=None,u=None):
    """ Linear mass matrix with small angle approximation
    q:  degrees of freedom at operating point, array-like: ['x(t)', 'y(t)', 'phi_x(t)', 'phi_y(t)', 'phi_z(t)', 'q_T1(t)', 'psi(t)']
    qd: dof velocities at operating point, array-like
    p:  parameters, dictionary with keys: ['L_T', 'tilt', 'v_yT1c', 'x_NR', 'z_NR']
    u:  inputs at operating point, dictionary with keys: []
           where each values is a constant!
    The columns of B correspond to:   [T_a(t)]\\ 
    """
    BB = np.zeros((7,1))
    BB[0,0] = -p['v_yT1c']*q[5]*sin(p['tilt'])-(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*q[3]+cos(p['tilt'])
    BB[1,0] = -(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))*q[4]+(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*q[2]
    BB[2,0] = -(-p['L_T']*(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))+p['x_NR']*sin(p['tilt'])+p['z_NR']*cos(p['tilt'])+(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*q[5])*q[4]
    BB[3,0] = -p['L_T']*(p['v_yT1c']*q[5]*sin(p['tilt'])-cos(p['tilt']))+p['x_NR']*sin(p['tilt'])+p['z_NR']*cos(p['tilt'])+(p['v_yT1c']*q[5]*cos(p['tilt'])+sin(p['tilt']))*q[5]
    BB[4,0] = 0
    BB[5,0] = p['v_yT1c']*p['x_NR']*sin(p['tilt'])+p['v_yT1c']*p['z_NR']*cos(p['tilt'])-p['v_yT1c']*q[5]*sin(p['tilt'])+cos(p['tilt'])
    BB[6,0] = 0
    return BB

