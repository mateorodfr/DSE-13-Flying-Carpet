

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
    """

    #Motor List
    #The object with index 0 is the currently selected one. All other indices are for comparison
    motor0 = ["Siemens SP260D",8,0.418/2,1500,204000,49,0.95,0.8] #Strong motor used on City airbus
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

    #function that returns the object properties in an array
    #can also print properties if given True as an argument
    def getParameters(self,printParameters=False):

        self.Parameters_motor = [self.Name_motor, self.N_motor,self.R_motor,self.Torque,self.P_max, self.M_motor, self.eff_motor, self.SF_rotational]
        if printParameters:
            print(f'The current motor parameters are:'
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
    """

    #Propeller list
    #The object with index 0 is the currently selected one. All other indices are for comparison
    #For the width and thickness of blade please enter the ratio in terms of Diameter i.e W/D & t/D

    propeller0 = ['NACA2412', 1.3, 6, 2.5, 0.1,0.01,660,0.9]
    propeller1 = ['NACA2412', 1.3, 6, 3.5, 0.1, 0.01, 660, 0.9]

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

    def getParameters(self, printParameters=False):

        self.Parameters_propeller = [self.airfoil_prop, self.CL_prop, self.N_prop,self.D_prop,self.W_prop,self.t_prop, self.rho_prop, self.eff_prop, self.S_prop,self.M_blades]
        if printParameters:
            print(f'The current propeller parameters are:'
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

    #Battery list
    #The object with index 0 is the currently selected one. All other indices are for comparison
    #The battery cost is a constant currently set at 100$/kWh

    cost_density = 0.1 #$/Wh

    battery0 = ['Panasonic NCA Si-C', 260, 683e3,0.9]

    def __init__(self,key=0):

        if key == 0:
            self.Name_battery = self.battery0[0] # Get source
            self.rhoE_battery = self.battery0[1] #[Wh/kg]
            self.rhoV_battery = self.battery0[2] #[Wh/m3]
            self.eff_battery = self.battery0[3] #get source this is an assumption
            self.rhoC_battery = self.rhoV_battery*self.cost_density # 100 dollars/kWh this is what tesla would like before 2020 so not an actual value, range is usual 150-125

    def getParameters(self, printParameters=False):

        self.Parameters_battery = [self.Name_battery,self.rhoE_battery,self.rhoV_battery,self.rhoC_battery, self.eff_battery]
        if printParameters:
            print(f'The current battery parameters are:'
            + f'\n\tName of Battery: {self.Parameters_battery[0]} [-]'
            + f'\n\tSpecific Energy of Battery: {self.Parameters_battery[1]} [Wh/kg]'
            + f'\n\tVolumetric Density of Battery: {self.Parameters_battery[2]} [Wh/m3]'
            + f'\n\tCost density of Battery: {self.Parameters_battery[3]} [$/kg]'
            + f'\n\tEfficiency of Battery: {self.Parameters_battery[4]} [-]'
            + f'\nThe variables in the output array appear in this order in SI units.'
            )
        return self.Parameters_battery

motor = BatteryParameters(0)
print(motor.getParameters(True))