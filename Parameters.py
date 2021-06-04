import numpy as np


'''
Introduction to Parameters
This script contains the individual components of the design
When calling this class in a module it creates an object based on the given parameters
The classes are called using a key (integer from 0 to n) this key corresponds to a certain configuration
By default when calling the class without a kehy it initializes key 0 which is the current selected componenent.


The class works as such, calling the class creates an object with its properties.
This allows to call the ibject and its variables in a script.
Furthermroe the function getParameters() outputs all the parameters in an array.
Optionally you can pass on a True argument to the fucntion to print the parameters

Check the class syntax to get familiar with the order of variables and how to initialze a configuration.
Different configurations can be input by creating an array inside the class which contains all relevenat information to build the object


Example syntax:

motor = MotorParameters(0) -> Create a motor object based on configuration 0 (default)
motorParams = motor.getParameters(True) -> Stores the motor parameters in an array and print an informative message
P_eng = motor.P_max -> assigns the max motor power to a variable.

Example Initialization Array for Motor:

motor0 = ["Siemens SP260D",8,0.418/4*1.15,1500/2,204000/2,49/2*1.15,0.95,0.8]
motor1 = ["Siemens SP260D",4,0.418/2,1500,204000,49,0.95,0.8]


UPDATED
if you create a concept object this will return all relevant subclasses.
This means that by passing a key between 0-3 you can generate each different concept
These will contain the relevant conponent parameters for each configuration

'''

class MotorParameters(object):

    """
    MOTORS
    These variables are all related to physical properties of the current motors
    Name_motor: Name of the motor
    N_motor: Number of engines and thus the propellers. If assumed a quad copter, the engines are rotated by
                180 deg and attached on top
    R_motor: The diameter of the motor is found on internet. Reference by Chloe
    Torque: The maximum continuous torque is 1500 Nm. When the critical loading configuration has been determined
            the Torque can be adapted, such that the power usage is half of the maximum needed.
    P_max: The maximum engine power in watts, found from literature online.
    M_motor: Motor mass, found on internet [kg]
    SF_rotational: Safety factor for moment of inertia of rotational part of the engine

    Initialization array: [Name_motor, N_motor, R_motor, Torque, P_max, M_motor, SF_rotational]
    """

    #Motor List
    #The object with index 0 is the currently selected one. All other indices are for comparison
    motor0 = ["Siemens SP260D",8,0.418/4*1.15,1500/2,204000/2,49/2*1.15,0.95,0.8] #Strong motor used on City airbus
    motor1 = ["Siemens SP260D",4,0.418/2,1500,204000,49,0.95,0.8]

    #Initialization if no key is given takes the default object '0'
    def __init__(self,key=0):
        if key == 0:
            self.Name_motor = self.motor0[0]
            self.N_motor = self.motor0[1]
            self.R_motor = self.motor0[2]
            self.Torque = self.motor0[3]
            self.P_max = self.motor0[4]
            self.M_motor = self.motor0[5]
            self.eff_motor = self.motor0[6]
            self.SF_rotational = self.motor0[7]

        elif key == 1:
            self.Name_motor = self.motor1[0]
            self.N_motor = self.motor1[1]
            self.R_motor = self.motor1[2]
            self.Torque = self.motor1[3]
            self.P_max = self.motor1[4]
            self.M_motor = self.motor1[5]
            self.eff_motor = self.motor1[6]
            self.SF_rotational = self.motor1[7]

    #function that returns the object properties in an array
    #can also print properties if given True as an argument
    def getParameters(self,printParameters=False):

        self.Parameters_motor = [self.Name_motor, self.N_motor,self.R_motor,self.Torque,self.P_max, self.M_motor, self.eff_motor, self.SF_rotational]
        if printParameters:
            print(
              f'*******************************************************************************'
            + f'\nThe current MOTOR parameters are:'
            + f'\n\tName: {self.Parameters_motor[0]} [-]'
            + f'\n\tNumber of Motors: {self.Parameters_motor[1]} [-]'
            + f'\n\tRadius of Motor: {self.Parameters_motor[2]} [m]'
            + f'\n\tMaximum Motor Torque: {self.Parameters_motor[3]} [Nm]'
            + f'\n\tMaximum Motor Power: {self.Parameters_motor[4]/1000} [kW]'
            + f'\n\tMass of Motor: {self.Parameters_motor[5]} [kg]'
            + f'\n\tEfficiency of Motor: {self.Parameters_motor[6]} [-]'
            + f'\n\tSafety Factor for Motor Moment of Inertia: {self.Parameters_motor[7]} [-]'
            + f'\nThe variables in the output array appear in this order in SI units.'
            )
        return self.Parameters_motor

class PropellerParameters(object):

    """
    PROPELLERS
    CL: Lift coefficient The cross section is assumed to be constant throughout the blade and the airfoil applied is NACA-2412
        CL determined from airfoiltools by taking the heighest Reynolds number
    N_blade: The blade count is assumed for now. No papers were found which derive the optimal number. Also,
             it is not clear whether the amount of blades reduces the efficiency of the lift generation properties
    D_blade: Various diameters of the propellers are inspected to determine the best fit
    W_blade: The width of the blade is assumed to be a tenth of the length. Could search for papers on propeller design
    t_blade: Thickness of the blade assumed to be 1% of length on average
    rho_blade: density of the propeller blades
    eff_prop: efficieny of propeller

    Initialization array: [Airfoil_name, CL, N_blade, D_blade, W/D , t/D, rho_blade, eff_prop]
    """

    #Propeller list
    #The object with index 0 is the currently selected one. All other indices are for comparison
    #For the width and thickness of blade please enter the ratio in terms of Diameter i.e W/D & t/D

    propeller0 = ['NACA2412', 1.3, 2, 2.5, 0.1,0.01,660,0.9]
    propeller1 = ['NACA2412', 1.3, 4, 3.5, 0.1, 0.01, 660, 0.9]

    #Initialization if no key is given takes the default object '0'
    def __init__(self,key=0):
        if key == 0:

            self.airfoil_prop = self.propeller0[0]
            self.CL_prop = self.propeller0[1]
            self.N_prop = self.propeller0[2]
            self.D_prop = self.propeller0[3]
            self.W_prop = self.propeller0[4] * self.D_prop
            self.t_prop = self.propeller0[5] * self.D_prop / 2
            self.rho_prop = self.propeller0[6]
            self.eff_prop = self.propeller0[7]
            self.S_prop = self.W_prop * self.D_prop * self.N_prop
            self.M_blades = self.S_prop * self.t_prop * self.rho_prop

        if key == 1:

            self.airfoil_prop = self.propeller1[0]
            self.CL_prop = self.propeller1[1]
            self.N_prop = self.propeller1[2]
            self.D_prop = self.propeller1[3]
            self.W_prop = self.propeller1[4] * self.D_prop
            self.t_prop = self.propeller1[5] * self.D_prop / 2
            self.rho_prop = self.propeller1[6]
            self.eff_prop = self.propeller1[7]
            self.S_prop = self.W_prop * self.D_prop * self.N_prop
            self.M_blades = self.S_prop * self.t_prop * self.rho_prop

    def getParameters(self, printParameters=False):

        self.Parameters_propeller = [self.airfoil_prop, self.CL_prop, self.N_prop,self.D_prop,self.W_prop,self.t_prop, self.rho_prop, self.eff_prop, self.S_prop,self.M_blades]
        if printParameters:
            print(
              f'*******************************************************************************'
            + f'\nThe current PROPELLER parameters are:'
            + f'\n\tAirfoil: {self.Parameters_propeller[0]} [-]'
            + f'\n\tLift Coefficient of Propeller: {self.Parameters_propeller[1]} [-]'
            + f'\n\tNumber of Propeller: {self.Parameters_propeller[2]} [-]'
            + f'\n\tDiameter of Propeller: {self.Parameters_propeller[3]} [m]'
            + f'\n\tAverage Width of Propeller: {self.Parameters_propeller[4]} [m]'
            + f'\n\tAverage Thickness of Propeller: {self.Parameters_propeller[5]} [m]'
            + f'\n\tDensity of Propeller: {self.Parameters_propeller[6]} [kg/m3]'
            + f'\n\tEfficiency of Propeller: {self.Parameters_propeller[7]} [-]'
            + f'\n\tTotal Surface area of Propeller: {self.Parameters_propeller[8]} [m2]'
            + f'\n\tTotal Mass of Propeller: {self.Parameters_propeller[9]} [kg]'
            + f'\nThe variables in the output array appear in this order in SI units.'
            )
        return self.Parameters_propeller

class BatteryParameters(object):

    """
    BATTERY
    Name_battery: The NACA identification number of the chosen airfoil
    rhoE_battery: Specific energy density of battery in [Wh/kg]
    rhoV_battery: Volumetric energy density of battery in [Wh/m3]
    eff_battery: Efficiency of the battery chosen
    rhoC_battery: Specific cost of batteries in [$/kg] this is calculated using the average cost_density

    Initialization array: [Name_battery, rhoE_battery, rhoV_battery,DOD, eff_battery]
    """

    #Battery constants
    cost_density = 0.15 #$/Wh

    #Battery list
    #The object with index 0 is the currently selected one. All other indices are for comparison
    #The battery cost is a constant currently set at 100$/kWh
    battery0 = ['Panasonic NCA Si-C', 300, 683000,969,0.7,0.9, 0.95, 0.15, 0.95, 25]


    def __init__(self,key=0):

        if key == 0:
            self.Name_battery = self.battery0[0] # Get source
            self.rhoE_battery = self.battery0[1] #[Wh/kg]
            self.rhoV_battery = self.battery0[2] #[Wh/m3]
            self.rhoP_battery = self.battery0[3] #W/kg
            self.dod_battery = self.battery0[4]
            self.eff_battery = self.battery0[5] #get source this is an assumption
            self.rhoC_battery = self.rhoE_battery*self.cost_density # 100 dollars/kWh this is what tesla would like before 2020 so not an actual value, range is usual 150-125
            self.eff_inverter = self.battery0[6]
            self.loss_factor = self.battery0[7]
            self.degradation = self.battery0[8]
            self.inverter_mass = self.battery0[9]

    def getParameters(self, printParameters=False):

        self.Parameters_battery = [self.Name_battery,self.rhoE_battery,self.rhoV_battery,self.rhoP_battery,self.dod_battery,self.rhoC_battery, self.eff_battery]
        if printParameters:
            print(

              f'*******************************************************************************'
            + f'\nThe current BATTERY parameters are:'
            + f'\n\tName of Battery: {self.Parameters_battery[0]} [-]'
            + f'\n\tSpecific Energy of Battery: {self.Parameters_battery[1]} [Wh/kg]'
            + f'\n\tVolumetric Density of Battery: {self.Parameters_battery[2]} [Wh/m3]'
            + f'\n\tPower Density of Battery: {self.Parameters_battery[3]} [W/kg]'
            + f'\n\tDepth of Discharge of Battery: {self.Parameters_battery[4]} [-]'
            + f'\n\tCost density of Battery: {self.Parameters_battery[5]} [$/kg]'
            + f'\n\tEfficiency of Battery: {self.Parameters_battery[6]} [-]'
            + f'\nThe variables in the output array appear in this order in SI units.'
            )
        return self.Parameters_battery

class CabinParameters(object):


    """
    Cabin initaliazation array
    L_person: Length of one person from back to tip of knees when seating
    W_person: Width of person from shoulder to shoulder when seating
    H_person: Height of person when seating
    t_wall: Thickness of the cabin walls
    nRows: number of rows when looking from front (for 6 seats with 2 seats abreast you have nRows = 3)
    nCols: number of seats abreast
    SF: Safety factor for dimensions

    Initialization Array: [L_person, W_person, H_person, t_wall, nRows, nCols, SF_dimensions]

    """
    #Cabin constants
    W_person = 0.75
    L_person = 0.5
    H_person = 1.5

    #Cabin list
    #The object with index 0 is the currently selected one. All other indices are for comparison
    cabin0 = [0.5,0.75,1.5,0.0005,3,2,1.1]

    def __init__(self,key=0):

        if key == 0:
            self.L_cabin = self.cabin0[-3]*self.cabin0[0]*self.cabin0[-1]
            self.W_cabin = self.cabin0[-2]*self.cabin0[1]*self.cabin0[-1]
            self.H_cabin = self.cabin0[2]*self.cabin0[-1]
            self.t_cabin = self.cabin0[3]
            self.S_cabin = [self.L_cabin*self.W_cabin,self.L_cabin*self.H_cabin,self.W_cabin*self.H_cabin] #LW, LH, WH
            self.V_cabin = self.L_cabin*self.W_cabin*self.H_cabin

    def getParameters(self,printParameters = False, getOutput = True):

        self.Parameters_cabin = [self.L_cabin,self.W_cabin,self.H_cabin,self.t_cabin,self.S_cabin, self.V_cabin]
        if printParameters:
            print(
              f'*******************************************************************************'
            + f'\nThe current CABIN parameters are:'
            + f'\n\tLength of Cabin: {np.round(self.Parameters_cabin[0],3)} [m]'
            + f'\n\tWidth of Cabin: {np.round(self.Parameters_cabin[1],3)} [m]'
            + f'\n\tHeight of Cabin: {np.round(self.Parameters_cabin[2],3)} [m]'
            + f'\n\tThickness of Cabin: {np.round(self.Parameters_cabin[3],3)} [m]'
            + f'\n\tSurface are of Faces of Cabin: {np.round(self.Parameters_cabin[4],3)} array[m2,m2,m2]'
            + f'\n\tVolume of Cabin: {np.round(self.Parameters_cabin[5],3)} [m3]'
            + f'\nThe variables in the output array appear in this order in SI units.'
            )
        return self.Parameters_cabin

class ConceptParameters(object):

    """
    CONCEPT
    Name_concept: The name of the concept
    Mtot_concept: Total concept mass
    Mpay_concept: Maximum payload mass for this concept
    Mbat_concept: Total concept battery mass
    Vtot_concept: Total volume of  concept
    Vbat_concept: Total battery volume of concept

    Initialization array: [Name_concept, Mtot_concept, Mpay_concept, Mbat_concept, Vtot_concept, Vbat_concept]
    """

    #Battery list
    #The object with index 0 is the currently selected one. All other indices are for comparison
    #The battery cost is a constant currently set at 100$/kWh
    concept0 = ['Pickup & Release', 1970 + 131 + 30, 600, 460 + 40, 12.3, 0.325, 51]
    concept1 = ['Pickup & Release', 1454, 600, 221.82, 12.3, 0.325]

    def __init__(self,key=0):

        if key == 0:
            self.Name_concept = self.concept0[0]
            self.Mtot_concept = self.concept0[1]
            self.Mpay_concept = self.concept0[2]
            self.Mbat_concept = self.concept0[3]
            self.Vtot_concept = self.concept0[4]
            self.Vbat_concept = self.concept0[5]
            self.M_LV_bat = self.concept0[6]


        if key == 1:
            self.Name_concept = self.concept1[0]
            self.Mtot_concept = self.concept1[1]
            self.Mpay_concept = self.concept1[2]
            self.Mbat_concept = self.concept1[3]
            self.Vtot_concept = self.concept1[4]
            self.Vbat_concept = self.concept1[5]

        self.motor = MotorParameters(key)
        self.propeller = PropellerParameters(key)
        self.battery = BatteryParameters(key)
        self.cabin = CabinParameters(key)
        self.physics = PhysicalParameters()

    def getParameters(self,printParameters = False,printComponents=False):

        self.Parameters_concept = [self.Name_concept,self.Mtot_concept,self.Mpay_concept,self.Mbat_concept,self.Vtot_concept, self.Vbat_concept]
        if printParameters:
            print(
              f'*******************************************************************************'
            + f'\nThe current CONCEPT parameters are:'
            + f'\n\tName of Concept: {self.Parameters_concept[0]} [-]'
            + f'\n\tMass Total of Concept: {self.Parameters_concept[1]} [kg]'
            + f'\n\tMass Payload of Concept: {self.Parameters_concept[2]} [kg]'
            + f'\n\tMass Battery of Concept: {self.Parameters_concept[3]} [kg]'
            + f'\n\tVolume Total of Concept: {self.Parameters_concept[4]} [m3]'
            + f'\n\tVolume Battery of Concept: {self.Parameters_concept[5]} [m3]'
            + f'\nThe variables in the output array appear in this order in SI units.'
            + f'\nThis concept contains a motor, propeller, battery and cabin object'
            )
        if printComponents:
            _ = self.motor.getParameters(True)
            _ = self.propeller.getParameters(True)
            _ = self.battery.getParameters(True)
            _ = self.cabin.getParameters(True)
        return self.Parameters_concept

    def MMOI(self):
        

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

        Mpayload = self.Mpay_concept
        mOEW = self.Mtot_concept-Mpayload- self.motor.N_motor*(self.motor.M_motor+self.propeller.M_blades)
        Dim_cabin = np.array([self.cabin.L_cabin,self.cabin.W_cabin,self.cabin.H_cabin])

        """Compute Moment of Inertias"""

        #Cabin moment of inertia (based on hollow cube mass moment of inertia, axis through cg)
        I_default = (5/18)*mOEW*Dim_cabin**2

        #Mass momnet of inertia around x axis
        I_xx_engines = (self.motor.N_motor*(self.motor.M_motor+self.propeller.M_blades))*(self.cabin.L_cabin/2+self.propeller.D_prop/2)**2 #Engine contribution
        rpayloadx = (self.cabin.L_person/2 + (0.1 * self.cabin.L_cabin)/2) #Payload zy distance
        I_xx_payload = Mpayload*rpayloadx**2 #Payload Contribution
        I_xx = I_default[0] + I_xx_engines+ I_xx_payload #Total x mass moment of inertia

        #Mass momnet of inertia around y axis
        I_yy_engines = (self.motor.N_motor*(self.motor.M_motor+self.propeller.M_blades))*(self.cabin.W_cabin/2+self.propeller.D_prop/2)**2 #Engine contribution
        rpayloady = (self.cabin.W_person + (self.cabin.W_cabin - 3 * self.cabin.W_person)/4) #Payload zx distance
        I_yy_payload = (2/3)*Mpayload*rpayloady**2 #Payload Contribution
        I_yy = I_default[1] + I_yy_engines + I_yy_payload #Total y mass moment of inertia

        #Mass moment of inertia around z axis
        r2enginez = ((self.cabin.L_cabin+self.propeller.D_prop)/2)**2 + ((self.cabin.W_cabin+self.propeller.D_prop)/2)**2 #Distance to the engine in xy (squared already)
        I_zz_engines = (self.motor.N_motor*(self.motor.M_motor+self.propeller.M_blades))*r2enginez #Engine Contribution
        r2payload1256 = rpayloadx**2 + rpayloady**2 #Distance to seats 1256 in xy (squared already)
        I_zz_payload = (2/3)*Mpayload*r2payload1256 + (1/3)*Mpayload*rpayloady**2 #Payload Contribution
        I_zz = I_default[2]+I_zz_engines+I_zz_payload #Total mass moment of inertia in z

        I = np.array([I_xx,I_yy,I_zz]) #Array to store mass moments of inertia

        return I

    def Thrust(self,n):
        P_eng = np.arange(0,self.motor.P_max+1,self.motor.P_max/n).astype(np.float32) #Single motor power range
        T_eng = ((np.pi/2)*self.propeller.D_prop**2*self.physics.rho0*(P_eng*self.propeller.eff_prop)**2)**(1/3)/self.motor.eff_motor
        return T_eng,P_eng


class PhysicalParameters(object):

    def __init__(self):
        self.rho0 = 1.225
        self.g = 9.80665

class ElectronicsParameters(object):
    Electronics0 = [12, 1.3, 4, 30*10**-6, 0.005, 10, 356, 0.3, 8, 300, 0.4, 2, 150, 0.3, 2, 200, 0.4, 5, 25, 0.1, 6, 0.95, 0.95, 0.9, 130, 4]

    # Initialization if no key is given takes the default object '0'
    def __init__(self, key=0):
        if key == 0:
            self.camera_power = self.Electronics0[0]
            self.camera_mass = self.Electronics0[1]
            self.camera_amount = self.Electronics0[2]
            self.T_sens_power = self.Electronics0[3]
            self.T_sens_mass = self.Electronics0[4]
            self.T_sens_amount = self.Electronics0[5]
            self.motor_controller_power = self.Electronics0[6]
            self.motor_controller_mass = self.Electronics0[7]
            self.motor_controller_amount = self.Electronics0[8]
            self.VCU_power = self.Electronics0[9]
            self.VCU_mass = self.Electronics0[10]
            self.VCU_amount = self.Electronics0[11]
            self.FC_power = self.Electronics0[12]
            self.FC_mass = self.Electronics0[13]
            self.FC_amount = self.Electronics0[14]
            self.AMS_power = self.Electronics0[15]
            self.AMS_mass = self.Electronics0[16]
            self.AMS_amount = self.Electronics0[17]
            self.SN_power = self.Electronics0[18]
            self.SN_mass = self.Electronics0[19]
            self.SN_amount = self.Electronics0[20]
            self.eff_buck = self.Electronics0[21]
            self.eff_boost = self.Electronics0[22]
            self.PF_electronics = self.Electronics0[23]
            self.pump_power = self.Electronics0[24]
            self.pump_amount = self.Electronics0[25]




