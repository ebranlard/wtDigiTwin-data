import numpy as np
from welib.yams.flexibility import GeneralizedMCK_PolyBeam, GMBeam 

def FAST2WisdemInputs(FST_file):
    """ 
    Generates WISDEM inputs to the frequency domain component (structural part) from a OpenFAST model
    This is an example, and it will be handled by Wisdem directly. 

    Feel free to rename the keys of the output, but adapt `WisdomInputs2StructureInputs` accordingly.

    INPUTS:
      - FST_file: path to OpenFAST input file
    OUTPUTS:
      - pW: Wisdem inputs to frequency domain component (structural part: Twr+RNA).

    """
    # --------------------------------------------------------------------------------}
    # ---  All this will be handled by WISDEM
    # --------------------------------------------------------------------------------{
    import os
    try:
        import welib.weio as weio
    except:
        import weio

    def R_y(t):
        return np.asarray( [[np.cos(t),0,np.sin(t)], [0,1,0], [-np.sin(t),0,np.cos(t)] ])
    # --- Reading main OpenFAST files
    ext=os.path.splitext(FST_file)[1]
    FST=weio.read(FST_file)
    rootdir = os.path.dirname(FST_file)
    EDfile = os.path.join(rootdir,FST['EDFile'].strip('"')).replace('\\','/')
    ED      = weio.read(EDfile)
    rootdir = os.path.dirname(EDfile)
    bldfile = os.path.join(rootdir,ED['BldFile(1)'].strip('"')).replace('\\','/')
    twrfile = os.path.join(rootdir,ED['TwrFile'].strip('"')).replace('\\','/')

    # --- Nacelle geometry
    theta_tilt_y = -ED['ShftTilt']*np.pi/180  # NOTE: tilt has wrong orientation in FAST
    R_NS = R_y(theta_tilt_y)
    r_NS_inN    = np.array([0             , 0, ED['Twr2Shft']]) # Shaft start in N
    r_SR_inS    = np.array([ED['OverHang'], 0, 0             ]) # Rotor center in S
    r_SGhub_inS = np.array([ED['HubCM']   , 0, 0             ]) + r_SR_inS # Hub G in S
    r_NR_inN    = r_NS_inN + R_NS.dot(r_SR_inS)                 # Rotor center in N
    r_NGnac_inN = np.array([ED['NacCMxn'],0,ED['NacCMzn']    ]) # Nacelle G in N
    r_RGhub_inS = - r_SR_inS + r_SGhub_inS

    # --- Bld
    bld = weio.read(bldfile).toDataFrame()
    z   = bld['BlFract_[-]']*(ED['TipRad']-ED['HubRad']) + ED['HubRad']
    m   = bld['BMassDen_[kg/m]']
    nSpan = len(z)
    nShapes=0
    PhiU = np.zeros((nShapes,3,nSpan)) # Shape
    s_G  = np.zeros((3,nSpan))
    s_G[2,:] = z
    jxxG= z*0 + m # NOTE: unknown
    MM, Gr, Ge, Oe, Oe6 = GMBeam(s_G, z, m, PhiU, jxxG=jxxG, bUseIW=True, main_axis='z', split_outputs=False, rot_terms=True)
    L_B = ED['TipRad']-ED['HubRad'] # Blade length
    M_B = MM[0,0]  # Blade mass
    theta_cone_y = ED['Precone(1)']*np.pi/180

    # --- "Hub" = Hub + Shaft + Gen
    M_hub       = ED['HubMass']
    z_HG        = r_SGhub_inS[2]
    IR_hub      = np.zeros((3,3))
    IR_hub[0,0] = ED['HubIner'] + ED['GenIner']*ED['GBRatio']**2 # Hub inertia at R in S
    #IG_hub = translateInertiaMatrix(IR_hub, M_hub, np.array([0,0,0]), r_RGhub_inS) # Hub inertia at G_H in S

    # --- Rotor
    nB        = ED['NumBl']
    M_blds    = nB * L_B

    # --- Nacelle alone
    M_nac       = ED['NacMass']
    I0_nac      = np.zeros((3,3))
    I0_nac[2,2] = ED['NacYIner'] # Inertia of nacelle at N in N
    #IG_nac = translateInertiaMatrixToCOG(I0_nac,M_nac, -r_NGnac_inN) # Inertia of nacelle at G_N in N 

    # --- RNA
    # RNA is rotor + hub + shaft + nacelle + generator as rigid body
    M_RNA = M_blds + M_hub + M_nac
    IG_RNA  = I0_nac # TODO this is wrong, but will be handled by WISDEM

    # --- Twr
    twr    = weio.read(twrfile)
    twrProp = twr.toDataFrame()

    # --------------------------------------------------------------------------------}
    # --- Setting Wisdem parameters
    # --------------------------------------------------------------------------------{
    pW=dict()
    # Twr
    pW['TowerSpan']     = twrProp['HtFract_[-]'].values*(ED['TowerHt']-ED['TowerBsHt']) # from 0 to tower length
    pW['TowerMassDens'] = twrProp['TMassDen_[kg/m]'].values
    pW['TowerFAStiff']  = twrProp['TwFAStif_[Nm^2]'].values
    pW['TowerSSStiff']  = twrProp['TwSSStif_[Nm^2]'].values
    pW['TowerDamp']     = np.array([twr['TwrFADmp(1)'], twr['TwrFADmp(2)'], twr['TwrSSDmp(1)'], twr['TwrSSDmp(2)']]) # structural damping ratio 
    pW['TowerCoeffs']   = np.array([[ twr['TwFAM1Sh(2)'], twr['TwFAM2Sh(2)'], twr['TwSSM1Sh(2)'], twr['TwSSM2Sh(2)']],
                                  [ twr['TwFAM1Sh(3)'], twr['TwFAM2Sh(3)'], twr['TwSSM1Sh(3)'], twr['TwSSM2Sh(3)']],
                                  [ twr['TwFAM1Sh(4)'], twr['TwFAM2Sh(4)'], twr['TwSSM1Sh(4)'], twr['TwSSM2Sh(4)']],
                                  [ twr['TwFAM1Sh(5)'], twr['TwFAM2Sh(5)'], twr['TwSSM1Sh(5)'], twr['TwSSM2Sh(5)']],
                                  [ twr['TwFAM1Sh(6)'], twr['TwFAM2Sh(6)'], twr['TwSSM1Sh(6)'], twr['TwSSM2Sh(6)']]])

    pW['TowerExp']  = np.arange(2,7) # Exponents for shape functions polynomial, OpenFAST typically use 2,3,4,5,6
    pW['TowerHt']   = ED['TowerHt']
    pW['TowerBsHt'] = ED['TowerBsHt']
    pW['Gravity']   = ED['Gravity']

    # Fnd
    #pW['PtfmCMzt']  = ED['PtfmCMzt']
    #pW['PtfmMass']  = ED['PtfmMass']
    #pW['PtfmRIner'] = ED['PtfmRIner']
    #pW['PtfmPIner'] = ED['PtfmPIner']
    #pW['PtfmYIner'] = ED['PtfmYIner']

    # RNA
    pW['ShftTilt']  = ED['ShftTilt']  
    pW['x_NR']      = r_NR_inN[0]     # x-coord from N to R in nac-coord
    pW['z_NR']      = r_NR_inN[2]     # z-coord from N to R in nac-coord
    pW['x_RNAG']    = r_SGhub_inS[0]  # x-coord from N to RNA_COG in nac-coord
    pW['z_RNAG']    = r_SGhub_inS[2]  # z-coord from N to RNA_COG in nac-coord
    pW['M_RNA']     = M_RNA           # Total mass of RNA
    pW['J_xx_RNA']  = IG_RNA[0,0]     # Inertia of RNA at RNA_COG in nac-coord # NOTE diagonal for now..
    pW['J_yy_RNA']  = IG_RNA[1,1]     # Inertia of RNA at RNA_COG in nac-coord
    pW['J_zz_RNA']  = IG_RNA[2,2]     # Inertia of RNA at RNA_COG in nac-coord
    return pW

def WisdemInputs2StructureInputs(pW):
    """ 
    Extract necessary information from Wisdom inputs and setup inputs required by linear structural model

    If possible do no change the keys used for the outputs `p`. 

    INPUTS: 
      pW : dictionary of wisdem inputs
    OUTPUTS
      p : dictionary for linear structural model
    
    """ 
    # --- Compute necessary tower flexibility parameters:
    # Generalized Mass, Stiffness, and Damping matrices, Shape integrals, shape function slopes
    s_span    = pW['TowerSpan']
    m         = pW['TowerMassDens']
    EIFA      = pW['TowerFAStiff']
    EISS      = pW['TowerSSStiff']
    coeffs    = pW['TowerCoeffs']
    exp       = pW['TowerExp']
    damp_zeta = pW['TowerDamp']/100
    gravity   = pW['Gravity']
    Mtop      = pW['M_RNA']
    coeffs    = coeffs[:,:1]  # NOTE: taking the first shape function only for now

    pTwr = GeneralizedMCK_PolyBeam(s_span, m, EIFA, EISS, coeffs, exp, damp_zeta, gravity=gravity, Mtop=Mtop, bAxialCorr=False, bStiffening=True)

    # --- Dict needed by structural script 
    p = dict()
    #p['z_FG']     = pW['PtfmCMzt']
    #p['M_F']      = pW['PtfmMass']
    #p['J_xx_F']   = pW['PtfmRIner']
    #p['J_yy_F']   = pW['PtfmPIner']
    #p['J_zz_F']   = pW['PtfmYIner']
    p['g']        = pW['Gravity']
    p['tilt']     =-pW['ShftTilt']
    p['x_NR']     = pW['x_NR']                    # x-coord from N to R in nac-coord
    p['z_NR']     = pW['z_NR']                    # z-coord from N to R in nac-coord
    p['x_RNAG']   = pW['x_RNAG']                  # x-coord from N to RNA_G in nac-coord
    p['z_RNAG']   = pW['z_RNAG']                  # z-coord from N to RNA_G in nac-coord
    p['M_RNA']    = pW['M_RNA']                   # Total mass of RNA
    p['J_xx_RNA'] = pW['J_xx_RNA']                # Inertia of RNA at RNA_G in nac-coord
    p['J_yy_RNA'] = pW['J_xx_RNA']                # Inertia of RNA at RNA_G in nac-coord
    p['J_zz_RNA'] = pW['J_xx_RNA']                # Inertia of RNA at RNA_G in nac-coord
    p['L_T']      = pW['TowerHt']-pW['TowerBsHt']
    p['z_OT']     = pW['TowerBsHt']               # distance from "Origin" (MSL) to tower base
    p['M_T']      = pTwr['MM'][0,0]
    p['z_TG']     = pTwr['s_OG'][2]
    p['J_xx_T']   = pTwr['J_G'][0,0]
    p['J_yy_T']   = pTwr['J_G'][1,1]
    p['J_zz_T']   = pTwr['J_G'][2,2]
    p['MM_T']     = pTwr['MM']
    p['Oe_T']     = pTwr['Oe6']
    p['Gr_T']     = pTwr['Gr']
    p['Ge_T']     = pTwr['Ge']
    p['v_yT1c']   = pTwr['alpha'][1,0] # NOTE: will need to be adapted for more modes
    p['DD_T']     = pTwr['DD']
    p['KK_T']     = pTwr['KK']
    return p



