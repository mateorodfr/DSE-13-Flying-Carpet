import numpy as np

""" Constants and mission parameters """

# Motors
N_motor = 8                 # Number of engines and thus the propellers. If assumed a quad copter, the engines are rotated by
                            # 180 deg and attached on top
R_motor = 0.418 / 2         # The diameter of the motor is found on internet. Reference by Chloe
Torque = 1500               # The maximum continuous torque is 1500 Nm. When the critical loading configuration has been determined
                            # the Torque can be adapted, such that the power usage is half of the maximum needed. P_max_motor = 204 kW
M_eng = 49                  # kg. Motor mass, found on internet
SF_Rotational = 0.8         # Rotational part of the engine

# Propellers
CL = 1.3                    # The cross section is assumed to be constant throughout the blade and the airfoil applied is NACA-2412
                            # CL determined from airfoiltools by taking the heighest Reynolds number
N_blade = 6                 # The blade count is assumed for now. No papers were found which derive the optimal number. Also,
                            # it is not clear whether the amount of blades reduces the efficiency of the lift generation properties
D_blade = 2.5    # Various diameters of the propellers are inspected to determine the best fit
W_blade = 0.1 * D_blade     # The width of the blade is assumed to be a tenth of the length. Could search for papers on propeller design
t_blade = 0.01 * D_blade / 2        # Thickness of the blade assumed to be 1% of length on average
rho_blade = 660                     # kg/m^3. Assumed to be be black walnut

# Mission
g = 9.81                    # kg/s^2
Maxpayload = 600            # 6 people, each around 100kg
M_person = 100              # kg
t_mission = 30/60           # TBD, a place holder as for now due to one cycle being of undetermined length. Expressed in hours
rho = 1.225                 # kg/m^3 - air density at sea level
visc = 1.46*10**-5              # Viscosity at sea level


# Battery
Bat_E_dense = 250           # Wh/kg A single high voltage battery assumed for now. A low voltage battery used for electronics
                            # However it would be less heavy, thus not included for now (As the power requirements are not known yet)
rho_battery = 400 * 10**3                       # Wh/m^3 Assumed for the time being as range is approximately 250 - 670

# Cabin
t_wall = 0.005
L_person = 0.5
W_person = 0.75
H_engine = 0.3
H_person = 1.5

""" Efficiencies of the propulsion and power subsystem """

Prop_eff = 0.85             # Propeller efficiency assumed for now
Motor_eff = 0.85            # Motor efficiency. An assumption
Bat_eff = 0.9               # Assumed battery efficiency


"""  Initial power and mass estimation of the craft  """
# Engine dimensions
Sb = W_blade * D_blade * N_blade                    # m^2 - lift area of all the blades (Only 2D area) for a single engine
M_blades = Sb * t_blade * rho_blade                 # Mass of all the blades in kg for a single engine
I_prop = 1/3 * M_blades * (D_blade*D_blade) / 4     # Mass moment of inertia for all the blades (A rectangle rotates around its base) of single motor
I_mot = 1/2 * M_eng * R_motor**2 * SF_Rotational    # Mass moment of inertia for an engine (A rotating disc around its centre)
I_tot = I_mot + I_prop                              # kg*m^2


# Engine usage power and propeller velocities
omega = Torque/I_tot                                # Angular velocity from angular momentum theory
V_tip = D_blade /2 * omega                         # Velocity of blades tips
P_eng = Torque * omega                              # Single engines power without efficiencies
print("P_eng for hovering = ", I_prop)

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
M_cabin = Mass_tot - Maxpayload
print("Total mass = ", Mass_tot)
print("Safety margin = ", T/Mass_tot)

"""   Initial area and volume sizing   """

SF_dimensions = 1.1                             # A assumed safety factor for the dimensions of the cabin to provide clearances
SF_engine_diemensions = 1.15

S_eng = np.pi * D_blade*D_blade / 4             # Area of an engine (Assuming a quadcopter design)
                                                # Value found on internet. Chloe has the reference
V_engine = S_eng * H_engine * SF_engine_diemensions   # Volume of a single engine (propeller + motor)
H_cabin = H_person * SF_dimensions              # Height of the cabin determined by seating configuration
L_cabin = L_person * 3 * SF_dimensions          # Determined from seating configuration
W_cabin = W_person * 2 * SF_dimensions          # Determined from seating configuration
S_cabin_bottom = W_cabin * L_cabin
S_cabin_side = H_cabin * L_cabin
S_cabin_front = W_cabin * H_cabin
V_cabin = H_cabin * W_cabin * L_cabin               # m^3
V_battery = M_bat * Bat_E_dense / rho_battery       # m^3
V_tot = V_battery + V_cabin + N_motor * V_engine    # m^3

print("Total Volume = ", V_tot)


"""   c.g computation and estimation  """

cg_init = np.array([L_cabin/2, W_cabin/2, H_cabin/2])                                       # C.g cabin
cg_seat_12 = np.array([(L_cabin - L_person/2 - (L_cabin-3*L_person)/4), W_cabin/2, H_cabin/2])    # C.g at longitudinal position
cg_seat_135 = np.array([L_cabin/2, W_person/2, H_cabin/2])
cg_tot_12 = (M_cabin * cg_init + 2 * M_person * cg_seat_12) / (M_cabin + 2 * M_person)
print(np.round(cg_tot_12, 3))
cg_tot_135 = (M_cabin * cg_init + 3 * M_person * cg_seat_135) / (M_cabin + 3 * M_person)
print(np.round(cg_tot_135, 3))
cg_range = (cg_tot_12 - cg_tot_135) * 2
print(cg_range)

Moment_seat12 = - N_motor * T / 2 * cg_range[0]
Moment_seat135 = - N_motor * T / 2 * cg_range[1]

print("Moment c.g longitudinal = ", Moment_seat12, "Nm")
print("Moment c.g lateral = ", Moment_seat135, "Nm")

T_diff_seat12 = Moment_seat12 / (cg_tot_12[0]+D_blade/2*SF_dimensions)
T_diff_seat135 = Moment_seat135 / (cg_tot_135[1]++D_blade/2*SF_dimensions)
print("T difference seat12 = ", T_diff_seat12)
print("T difference seat135 = ", T_diff_seat135)

Increment_seat12 = cg_range[0] / 2 / L_cabin
Increment_seat135 = cg_range[1] / 2 / W_cabin
print("The approximate increase in c.g due to passengers laterally", Increment_seat12 * 100, "%")
print("The approximate increase in c.g due to passengers longitudinally", Increment_seat135 * 100, "%")


"""   One engine pair inoperative   """
# The assumption is that the engines are paired having a rotation of 180deg between them. Thus, the net torque generated
# in case of one pair of engines inoperative, is still equal to zero. Thus only roll of the craft is present due to
# an upwards force. Thus spin of the motor would have to be altered continuously in a sinusoidal manner such that the
# roll angles are limited.


# Take a look at TU Delft method of drone stabilisation
# Discuss verbally the concept of one engine inoperative, as that could be counteracted by shutting down a diagonal engine
# Thus, can be discussed qualitatively

I_yy = 1/3 * (L_cabin**2 * W_cabin * t_wall)    # The wall thickness is assumed. The cross section is still TBD
I_xx = 1/3 * (W_cabin**2 * L_cabin * t_wall)
Operation_percentage = np.arange(0.2, 1, 0.1)   # Assume the engine is still providing thrust just after the failure,
                                                # however, the percentage of thrust provision is f(t). If there is a
                                                # large scale debri impact, the operational_perc = 0, thus, catastrophic event.
                                                #

My = (1 - Operation_percentage) * T * 2 * (L_cabin/2 + D_blade/2)
Mx = T * 2 * (W_cabin/2 + D_blade/2)

t_response = 0.0005                             # Unfeasible due to gyroscope limitations. Arvis has source for that

alpha_y_init = My/I_xx                          # Rotation can be assumed to be accelerated in no wind conditions.
alpha_x_init = Mx/I_yy                          # As the airflow is assumed to be zero for the critical conditions, the
                                                # change in angular acceleration can be assumed constant due the lack of
                                                # external forces.

theta_reached_x = alpha_x_init * t_response**2 * 57.3 / 2
theta_reached_y = alpha_y_init * t_response**2 * 57.3 / 2

print(alpha_x_init, theta_reached_y)


"""    Gust load estimation    """




