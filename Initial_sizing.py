import numpy as np

""" Constants and mission parameters """

# Motors
N_motor = 8                 # Number of engines and thus the propellers. If assumed a quad copter, the engines are rotated by
                            # 180 deg and attached on top
R_motor = 0.418 / 2         # The diameter of the motor is found on internet. Reference by Chloe
Torque = 1000               # The maximum continuous torque is 1500 Nm. When the critical loading configuration has been determined
                            # the Torque can be adapted, such that the power usage is half of the maximum needed. P_max_motor = 204 kW
M_eng = 49                  # kg. Motor mass, found on internet

# Propellers
CL = 0.9                    # The cross section is assumed to be constant throughout the blade and the airfoil applied is NACA-2412
                            # CL determined from airfoiltools by taking the heighest Reynolds number
N_blade = 6                 # The blade count is assumed for now. No papers were found which derive the optimal number. Also,
                            # it is not clear whether the amount of blades reduces the efficiency of the lift generation properties
D_blade = np.arange(2.5, 3, 0.1)    # Various diameters of the propellers are inspected to determine the best fit
W_blade = 0.1 * D_blade     # The width of the blade is assumed to be a tenth of the length. Could search for papers on propeller design
t_blade = 0.01 * D_blade / 2        # Thickness of the blade assumed to be 1% of length on average
rho_blade = 660                     # kg/m^3. Assumed to be be black walnut

# Mission
g = 9.81                    # kg/s^2
Maxpayload = 600            # 6 people, each around 100kg
t_mission = 45/60           # TBD, a place holder as for now due to one cycle being of undetermined length. Expressed in hours
rho = 1.225                 # kg/m^3 - air density at sea level
visc = 1.46*10**-5              # Viscosity at sea level


# Battery
Bat_E_dense = 250           # Wh/kg A single high voltage battery assumed for now. A low voltage battery used for electronics
                            # However it would be less heavy, thus not included for now (As the power requirements are not known yet)
rho_battery = 400 * 10**3                       # Wh/m^3 Assumed for the time being as range is approximately 250 - 670

""" Efficiencies of the propulsion and power subsystem """

Prop_eff = 0.85             # Propeller efficiency assumed for now
Motor_eff = 0.85            # Motor efficiency. An assumption
Bat_eff = 0.9               # Assumed battery efficiency


"""  Initial power and mass estimation of the craft  """
# Engine dimensions
Sb = W_blade * D_blade * N_blade                    # m^2 - lift area of all the blades (Only 2D area) for a single engine
M_blades = Sb * t_blade * rho_blade                 # Mass of all the blades in kg for a single engine
I_prop = 1/3 * M_blades * (D_blade*D_blade) / 4     # Mass moment of inertia for all the blades (A rectangle rotates around its base) of single motor
I_mot = 1/2 * M_eng * R_motor**2                    # Mass moment of inertia for an engine (A rotating disc around its centre)
I_tot = I_mot + I_prop                              # kg*m^2


# Engine usage power and propeller velocities
omega = Torque/I_tot                                # Angular velocity from angular momentum theory
V_tip = D_blade / 2 * omega                         # Velocity of blades tips
P_eng = Torque * omega                              # Single engines power without efficiencies
print("P_eng for hovering = ", P_eng)

# Reynolds number
V = V_tip / np.sqrt(2)              # It is assumed that the average velocity of the engine blades is the RMS value of the tip velocity. Needs to be changed possibly
Re = rho * V * W_blade / visc       # Reynolds number of the blades at sea level


# Thrust per engine
T = V_tip**2 / 6 * rho * CL * Sb * Prop_eff         # Thrust of a single engine. Only propeller efficiency taken into account as it corresponds to generation of mechanical energy
print("Mass to be carried by engines = ", T * N_motor / g)

# Power required and the mass estimations
P_req = P_eng * N_motor / Motor_eff / Prop_eff                          # Power required from the batteries for operative engines
M_bat = P_req * t_mission / Bat_E_dense / Bat_eff                       # Mass of the batteries
Mass_tot = 1.2 * (M_bat + Maxpayload + N_motor * (M_eng + M_blades))    # The structures mass is assumed to be 1.2 of the total mass, thus the introduction of the term

print("Total mass = ", Mass_tot)
print("Safety margin = ", T/Mass_tot)

"""   Initial area and volume sizing   """

SF_dimensions = 1.1                     # A assumed safety factor for the dimensions of the cabin to provide clearances
SF_engine_diemensions = 1.15

S_eng = np.pi * D_blade*D_blade / 4             # Area of an engine (Assuming a quadcopter design)
H_engine = 0.3                                  # Value found on internet. Chloe has the reference
V_engine = S_eng * H_engine * SF_dimensions     # Volume of a single engine (propeller + motor)
H_cabin = 1.5 * SF_dimensions                   # Height of the cabin determined by seating configuration
W_cabin = 0.75 * 3 * SF_dimensions              # Determined from seating configuration
L_cabin = 0.5 * 2 * SF_dimensions               # Determined from seating configuration
S_cabin_bottom = W_cabin * L_cabin
S_cabin_side = H_cabin * L_cabin
S_cabin_front = W_cabin * H_cabin
V_cabin = H_cabin * W_cabin * L_cabin               # m^3
V_battery = M_bat * Bat_E_dense / rho_battery       # m^3
V_tot = V_battery + V_cabin + N_motor * V_engine    # m^3

print("Total Volume = ", V_tot)


