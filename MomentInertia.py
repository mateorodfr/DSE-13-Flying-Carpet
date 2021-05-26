import numpy as np
import Parameters as pm


def MMOI(concept):


    """
    Mass Moment of Inertia function
    This fucntion takes a certain motor, propeller, cabin and concept objects
    This functions returns the MMOI in x,y,z plane in an array
    This fucntion shoudl work for any motor, propeller, cabin and concept combination

    Example syntax:
    motor = pm.MotorParameters(0)
    propeller = pm.PropellerParameters(0)
    cabin = pm.CabinParameters(0)
    concept = pm.ConceptParameters(0)
    I = MMOI(motor,propeller,cabin,concept)
    """

    """Initialize Parameters"""

    Mpayload = concept.Mpay_concept
    mOEW = concept.Mtot_concept-Mpayload- concept.motor.N_motor*(concept.motor.M_motor+concept.propeller.M_blades)
    Dim_cabin = np.array([concept.cabin.L_cabin,concept.cabin.W_cabin,concept.cabin.H_cabin])

    """Compute Moment of Inertias"""

    #Cabin moment of inertia (based on hollow cube mass moment of inertia, axis through cg)
    I_default = (5/18)*mOEW*Dim_cabin**2

    #Mass momnet of inertia around x axis
    I_xx_engines = (concept.motor.N_motor*(concept.motor.M_motor+concept.propeller.M_blades))*(concept.cabin.L_cabin/2+concept.propeller.D_prop/2)**2 #Engine contribution
    rpayloadx = (concept.cabin.L_person/2 + (0.1 * concept.cabin.L_cabin)/2) #Payload zy distance
    I_xx_payload = Mpayload*rpayloadx**2 #Payload Contribution
    I_xx = I_default[0] + I_xx_engines+ I_xx_payload #Total x mass moment of inertia

    #Mass momnet of inertia around y axis
    I_yy_engines = (concept.motor.N_motor*(concept.motor.M_motor+concept.propeller.M_blades))*(concept.cabin.W_cabin/2+concept.propeller.D_prop/2)**2 #Engine contribution
    rpayloady = (concept.cabin.W_person + (concept.cabin.W_cabin - 3 * concept.cabin.W_person)/4) #Payload zx distance
    I_yy_payload = (2/3)*Mpayload*rpayloady**2 #Payload Contribution
    I_yy = I_default[1] + I_yy_engines + I_yy_payload #Total y mass moment of inertia

    #Mass moment of inertia around z axis
    r2enginez = ((concept.cabin.L_cabin+concept.propeller.D_prop)/2)**2 + ((concept.cabin.W_cabin+concept.propeller.D_prop)/2)**2 #Distance to the engine in xy (squared already)
    I_zz_engines = (concept.motor.N_motor*(concept.motor.M_motor+concept.propeller.M_blades))*r2enginez #Engine Contribution
    r2payload1256 = rpayloadx**2 + rpayloady**2 #Distance to seats 1256 in xy (squared already)
    I_zz_payload = (2/3)*Mpayload*r2payload1256 + (1/3)*Mpayload*rpayloady**2 #Payload Contribution
    I_zz = I_default[2]+I_zz_engines+I_zz_payload #Total mass moment of inertia in z

    I = np.array([I_xx,I_yy,I_zz]) #Array to store mass moments of inertia

    return I


