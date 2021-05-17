

class MotorParameters(object):

    def __init__(self):

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


        self.Name_motor = "Siemens SP260D"
        self.N_motor = 8
        self.R_motor = 0.418 / 2
        self.Torque = 1500
        self.P_max = 204000
        self.M_motor = 49
        self.eff_motor = 0.95
        self.SF_rotational = 0.8

    def getMotorParameters(self,printParameters=False):

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
            + f'\nThe variables in the output array appear in this order.'
            )
        return self.Parameters_motor

motor = MotorParameters()
print(motor.getMotorParameters(True))

def PropellerParameters(object):

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

    self.CL = 1.3
    self.N_blade = 6
    self.D_blade = 2.5
    self.W_blade = 0.1 * self.D_blade
    self.t_blade = 0.01 * self.D_blade / 2
    self.rho_blade = 660