import numpy as np
import Parameters as pm

""" Constants and mission parameters """

# Motors
motor = pm.MotorParameters(0)
N_motor = motor.N_motor                 # Number of engines and thus the propellers. If assumed a quad copter, the engines are rotated by
                            # 180 deg and attached on top
R_motor = motor.R_motor        # The diameter of the motor is found on internet. Reference by Chloe
Torque = motor.Torque               # The maximum continuous torque is 1500 Nm. When the critical loading configuration has been determined
                            # the Torque can be adapted, such that the power usage is half of the maximum needed. P_max_motor = 204 kW
M_eng = motor.M_motor                  # kg. Motor mass, found on internet
SF_Rotational = motor.SF_rotational         # Rotational part of the engine
Motor_eff = motor.eff_motor

# Propellers
prop = pm.PropellerParameters(1)
CL = prop.CL_prop                    # The cross section is assumed to be constant throughout the blade and the airfoil applied is NACA-2412
                            # CL determined from airfoiltools by taking the heighest Reynolds number
N_blade = prop.N_prop                 # The blade count is assumed for now. No papers were found which derive the optimal number. Also,
                            # it is not clear whether the amount of blades reduces the efficiency of the lift generation properties
D_blade = prop.D_prop               # Various diameters of the propellers are inspected to determine the best fit
W_blade = prop.W_prop     # The width of the blade is assumed to be a tenth of the length. Could search for papers on propeller design
t_blade = prop.t_prop        # Thickness of the blade assumed to be 1% of length on average
rho_blade = prop.rho_prop                     # kg/m^3. Assumed to be be black walnut
Sb = prop.S_prop                    # m^2 - lift area of all the blades (Only 2D area) for a single engine
M_blades = prop.M_blades            # Mass of all the blades in kg for a single engine
Prop_eff = prop.eff_prop             # Propeller efficiency assumed for now


#Battery
battery = pm.BatteryParameters(0)
Bat_E_dense = battery.rhoE_battery           # Wh/kg A single high voltage battery assumed for now. A low voltage battery used for electronics
                            # However it would be less heavy, thus not included for now (As the power requirements are not known yet)
rho_battery = battery.rhoV_battery                       # Wh/m^3 Assumed for the time being as range is approximately 250 - 670
Bat_eff = battery.eff_battery
# Mission
g = 9.81                    # kg/s^2
Maxpayload = 600            # 6 people, each around 100kg
M_person = 100              # kg
#t_mission = 30/60          # TBD, a place holder as for now due to one cycle being of undetermined length. Expressed in hours
rho = 1.225                 # kg/m^3 - air density at sea level
visc = 1.46*10**-5          # Viscosity at sea level
t_ascend = 50/3600
t_hover = 10.5/60
t_desc = 86/3600
t_switch = 1/60 * 1                            # Multiplication by 0 means that external power line is present. Durring attachment some power is going to be taken

t_desc_acc = 5                                # s. The acceleration time. It can be assumed to be anything to have a decent complete ascend/descend time. Also, not too high loads
t_asc_acc = 10                                # s.  The acceleration time. It can be assumed to be anything to have a decent complete ascend/descend time. Also, not too high loads
h = 400

# Cabin
t_wall = 0.005
L_person = 0.5
W_person = 0.75
H_engine = 0.3
H_person = 1.5

""" Efficiencies of the propulsion and power subsystem """


            # Motor efficiency. An assumption



"""  Initial power and mass estimation of the craft  """
# Engine usage power and propeller velocities
#max_omega = Torque/I_tot                                # Angular velocity from angular momentum theory
P_eng = np.arange(10_000, 205_001, 5_000, dtype = np.int64)
#print("Initial ascend power = ", P_eng[12])
omega = P_eng/Torque * Motor_eff
V_tip = D_blade / 2 * omega                         # Velocity of blades tips
#print("V_tip = ", V_tip)
#P_eng = Torque * omega                              # Single engines power without efficiencies
#print("omega = ", omega)

# Reynolds number
V = V_tip / np.sqrt(2)              # It is assumed that the average velocity of the engine blades is the RMS value of the tip velocity. Needs to be changed possibly
Re = rho * V * W_blade / visc       # Reynolds number of the blades at sea level


# Thrust per engine
T = V_tip**2 / 6 * rho * CL * Sb * Prop_eff         # Thrust of a single engine. Only propeller efficiency taken into account as it corresponds to generation of mechanical energy

T2 = (np.pi/2 * D_blade**2 * rho * P_eng*P_eng * Motor_eff**2 * Prop_eff**2) ** (1/3)
#print("Thrust difference = ", T/T2)
#print("Mass possible to be carried by engines = ", (T2 / g))

# Ascend
t_asc_deac = t_asc_acc
h_asc_acc = (1.1 - 1) * g * t_asc_acc**2 / 2     # The deacceleration is assumed to be the same as the acceleration for those intervals
h_asc_deac = - (0.9 - 1) * g * t_asc_deac**2 / 2
h_left_asc = h - h_asc_acc - h_asc_deac
t_left_asc = h_left_asc/((1.1 - 1) * g * t_asc_acc)
t_asc_act = t_left_asc + 2 * t_asc_deac

# Descend
t_desc_deac = t_desc_acc
h_desc_deac = (1.1 - 1) * g * t_desc_acc**2 / 2     # The deacceleration is assumed to be the same as the acceleration for those intervals
h_desc_acc = - (0.9 - 1) * g * t_desc_acc**2 / 2
h_left = h - h_desc_acc - h_desc_deac
t_left = h_left/((- 0.9 + 1) * g * t_desc_acc)
t_desc_act = t_left + 2 * t_desc_deac

# Power required and the mass estimations
M_bat_asc_1 = N_motor * ((P_eng[8] + P_eng[11]) * t_asc_acc + P_eng[9] * t_left_asc) / Bat_E_dense /3600 / Bat_eff / Motor_eff / Prop_eff
M_bat_desc_1 = N_motor * ((P_eng[11] + P_eng[8]) * t_desc_acc + P_eng[9] * t_left) / Bat_E_dense /3600 / Bat_eff / Motor_eff / Prop_eff
M_bat_hov_1 = N_motor * (t_hover * P_eng[9]) / Bat_E_dense / Bat_eff / Motor_eff / Prop_eff
M_bat_switch_1 = N_motor * (t_switch * P_eng[9]) / Bat_E_dense / Bat_eff / Motor_eff / Prop_eff


Mass_tot1 = 1.2 * (M_bat_desc_1 + M_bat_asc_1 + M_bat_hov_1 + M_bat_switch_1 + Maxpayload + N_motor * (M_eng + M_blades))
M_cabin = Mass_tot1 - Maxpayload

#E_desc_act = M_bat_desc_1 * Bat_E_dense
#E_asc_act = M_bat_asc_1 * Bat_E_dense
#E_asc = 1.2 * (M_bat_desc_1 + M_bat_asc_1 + M_bat_hov_1 + M_bat_switch_1 + Maxpayload/6 + N_motor * (M_eng + M_blades)) * g * h / Prop_eff / Motor_eff / 3600
#E_desc = 1.2 * (M_bat_desc_1 + M_bat_asc_1 + M_bat_hov_1 + M_bat_switch_1 + Maxpayload + N_motor * (M_eng + M_blades)) * g * h / Prop_eff / Motor_eff / 3600

SM = T2*N_motor/(Mass_tot1*g)                      # Assumed deacceleration
print("Safety margin = ", T2*N_motor/(Mass_tot1*g))

print("Approximate time of descend (actual) ", t_desc_act)     # The acceleration of descent is assumed to be the same as for ascend only reversed.
                                    # However, the thrust durations are different to have 2 different times

print("Approximate time of ascend (actual) ", t_asc_act)

#t_descend_act = np.sqrt(-4 / ((SM-1) * g) * h)        # Multiplied by 2 just to have acceleration only at half the time and then deacceleration
#print("Approximate time of descend (actual) ", t_descend_act)


"""   Initial area and volume sizing   """

M_bat = M_bat_asc_1 + M_bat_hov_1 + M_bat_desc_1 + M_bat_switch_1

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

# print("Total Volume = ", V_tot)

print("Iterated total mass = ", Mass_tot1)

#print("Mass of the batteries needed for ascend", M_bat_asc_1)

"""   c.g computation and estimation  """

cg_init = np.array([L_cabin/2, W_cabin/2, H_cabin/2])                                       # C.g cabin
cg_seat_12 = np.array([(L_cabin - L_person/2 - (L_cabin-3*L_person)/4), W_cabin/2, H_cabin/2])    # C.g at longitudinal position
cg_seat_135 = np.array([L_cabin/2, W_person/2, H_cabin/2])
cg_tot_12 = (M_cabin * cg_init + 2 * M_person * cg_seat_12) / (M_cabin + 2 * M_person)
#print(np.round(cg_tot_12, 3))
cg_tot_135 = (M_cabin * cg_init + 3 * M_person * cg_seat_135) / (M_cabin + 3 * M_person)
#print(np.round(cg_tot_135, 3))
cg_range = (cg_tot_12 - cg_tot_135) * 2
#print(cg_range)

Moment_seat12 = - N_motor * T2 / 2 * cg_range[0]
Moment_seat135 = - N_motor * T2 / 2 * cg_range[1]

#print("Moment c.g longitudinal = ", Moment_seat12, "Nm")
#print("Moment c.g lateral = ", Moment_seat135, "Nm")

T_diff_seat12 = Moment_seat12 / (cg_tot_12[0]+D_blade/2*SF_dimensions)
T_diff_seat135 = Moment_seat135 / (cg_tot_135[1]++D_blade/2*SF_dimensions)
#print("T difference seat12 = ", T_diff_seat12)
#print("T difference seat135 = ", T_diff_seat135)

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

I_yy = 1/12 * ((W_cabin**2 + H_cabin**2) * (Mass_tot1 - N_motor * (M_eng + M_blades) - Maxpayload)) + N_motor * (M_eng + M_blades) * (D_blade/2 + W_cabin/2)**2 + 4 * Maxpayload / 6 * (W_person + (W_cabin - 3 * W_person)/4)**2 # The wall thickness is assumed. The cross section is still TBD
I_xx = 1/12 * ((L_cabin**2 + H_cabin**2) * (Mass_tot1 - N_motor * (M_eng + M_blades) - Maxpayload)) + N_motor * (M_eng + M_blades) * (D_blade/2 + L_cabin/2)**2 + Maxpayload * (L_person + (0.1 * L_cabin)/2)**2
Operation_percentage = np.array([1, 0.8, 0.6, 0.4, 0.2, 0.0])   # Assume the engine is still providing thrust just after the failure,
                                                # however, the percentage of thrust provision is f(t). If there is a
                                                # large scale debri impact, the operational_perc = 0, thus, catastrophic event.

My_init = (1 - Operation_percentage[-1]) * T2[3] * 2 * (L_cabin/2 + D_blade/2)
Mx_init = (1 - Operation_percentage[-1]) * T2[3] * 2 * (W_cabin/2 + D_blade/2)

t_response = 0.25                                   # Unfeasible due to gyroscope limitations.

alpha_y_init = My_init/I_xx                         # Rotation can be assumed to be accelerated in no wind conditions.
alpha_x_init = Mx_init/I_yy                         # As the airflow is assumed to be zero for the critical conditions, the
                                                    # change in angular acceleration can be assumed constant due the lack of
                                                    # external forces.
omega_x_init = alpha_x_init * t_response
omega_y_init = alpha_y_init * t_response

theta_reached_x = alpha_x_init * t_response**2 / 2
theta_reached_y = alpha_y_init * t_response**2 / 2


t_stop_engine = 0.5             # Engines are assumed to have the same velocity now. The functioning arm is now stopped and will be reversed to max thrust almost instantaneously
                                # A time in which an engine of this size can be stopped and reversed has to be determined!
theta_addition_x = theta_reached_x + omega_x_init * t_stop_engine
theta_addition_y = theta_reached_y + omega_y_init * t_stop_engine

t_reverse = 0.3             # Matter of choice

My_reverse = - (1 + Operation_percentage[-1]) * T2[7] / np.sqrt(2) * 2 * (L_cabin/2 + D_blade/2)
Mx_reverse = - (1 + Operation_percentage[-1]) * T2[7] / np.sqrt(2) * 2 * (W_cabin/2 + D_blade/2)

alpha_y_reverse = My_reverse/I_xx
alpha_x_reverse = Mx_reverse/I_yy

theta_zero_vec_x = theta_addition_x - omega_x_init**2 / alpha_x_reverse + alpha_x_reverse * (- omega_x_init / alpha_x_reverse)**2 / 2
theta_zero_vec_y = theta_addition_y - omega_y_init**2 / alpha_y_reverse + alpha_y_reverse * (- omega_y_init / alpha_y_reverse)**2 / 2

a_1_x_E = (N_motor * (3 + Operation_percentage[2]) / 4 * np.cos(np.pi / 2 - theta_reached_x) * (T2[3] + T2[4])/2) / Mass_tot1
a_1_y_E = (N_motor * (3 + Operation_percentage[2]) / 4 * np.cos(np.pi / 2 - theta_reached_y) * (T2[3] + T2[4])/2) / Mass_tot1
a_1_z_E = (- N_motor * (3 + Operation_percentage[2]) / 4 * np.cos(max(theta_reached_x, theta_reached_y) / np.sqrt(2)) * (T2[3] + T2[4])/2 + Mass_tot1 * g) / Mass_tot1

S_1_x_E = a_1_x_E / 2 * t_stop_engine**2
S_1_y_E = a_1_y_E / 2 * t_stop_engine**2
S_1_z_E = a_1_z_E / 2 * t_stop_engine**2

a_2_x_E = (N_motor * (2 + Operation_percentage[-1]) / 4 * np.cos(np.pi / 2 - theta_zero_vec_x / np.sqrt(2)) * T2[7] - N_motor / 4 * np.cos(np.pi / 2 - theta_zero_vec_x / np.sqrt(2)) * T2[7]) / Mass_tot1
a_2_y_E = (N_motor * (2 + Operation_percentage[-1]) / 4 * np.cos(np.pi / 2 - theta_zero_vec_y / np.sqrt(2)) * T2[7] - N_motor / 4 * np.cos(np.pi / 2 - theta_zero_vec_y / np.sqrt(2)) * T2[7]) / Mass_tot1
a_2_z_E = (- N_motor * (2 + Operation_percentage[-1]) / 4 * np.cos(max(theta_reached_x, theta_reached_y) / np.sqrt(2)) * T2[7] + N_motor / 4 * np.cos(max(theta_reached_x, theta_reached_y) / np.sqrt(2)) * T2[7] + Mass_tot1 * g) / Mass_tot1


theta_reverse_x = theta_addition_x + omega_x_init * t_reverse + alpha_x_reverse * t_reverse**2 / 2
theta_reverse_y = theta_addition_y + omega_y_init * t_reverse + alpha_y_reverse * t_reverse**2 / 2

V_1_x_E = a_1_x_E * t_stop_engine
V_1_y_E = a_1_y_E * t_stop_engine
V_1_z_E = a_1_z_E * t_stop_engine

S_x_E = S_1_x_E + V_1_x_E * t_reverse + a_2_x_E / 2 * t_reverse**2
S_y_E = S_1_y_E + V_1_y_E * t_reverse + a_2_y_E / 2 * t_reverse**2
S_z_E = S_1_z_E + V_1_z_E * t_reverse + a_2_z_E / 2 * t_reverse**2

#x_max = S_1_x_E - V_1_x_E**2 / a_2_x_E + (V_1_x_E / a_2_x_E)**2 / 2 * a_2_x_E
#y_max = S_1_y_E - V_1_y_E**2 / a_2_y_E + (V_1_y_E / a_2_y_E)**2 / 2 * a_2_y_E

V_x_E = V_1_x_E + a_2_x_E * t_reverse
V_y_E = V_1_y_E + a_2_y_E * t_reverse
V_z_E = V_1_z_E + a_2_z_E * t_reverse

#theta_tot_x = theta_addition_x + theta_reverse_x
#theta_tot_y = theta_addition_y + theta_reverse_y

omega_tot_x = omega_x_init + alpha_x_reverse * t_reverse
omega_tot_y = omega_y_init + alpha_y_reverse * t_reverse

t_stop_engine_2 = 0.5
theta_stop_x = theta_reverse_x + omega_tot_x * t_stop_engine_2      # Add the acc and velocity at this point
theta_stop_y = theta_reverse_y + omega_tot_y * t_stop_engine_2
t_stabilise_x = - 2 * theta_stop_x / omega_tot_x
t_stabilise_y = - 2 * theta_stop_y / omega_tot_y
alpha_stabilise_x = - omega_tot_x / t_stabilise_x
alpha_stabilise_y = - omega_tot_y / t_stabilise_y

M_x_needed = alpha_stabilise_x * I_yy
M_y_needed = alpha_stabilise_y * I_xx

T_needed_x = M_x_needed / 2 / (L_cabin/2 + D_blade/2)
T_needed_y = M_y_needed / 2 / (W_cabin/2 + D_blade/2)

a_3_x_E = (N_motor * (2 + Operation_percentage[-1]) / 4 * np.cos(np.pi / 2 - theta_zero_vec_x / np.sqrt(2)) * T2[7] + N_motor / 4 * np.cos(np.pi / 2 - theta_zero_vec_x / np.sqrt(2)) * T_needed_x) / Mass_tot1
a_3_y_E = (N_motor * (2 + Operation_percentage[-1]) / 4 * np.cos(np.pi / 2 - theta_zero_vec_y / np.sqrt(2)) * T2[7] + N_motor / 4 * np.cos(np.pi / 2 - theta_zero_vec_y / np.sqrt(2)) * T_needed_y) / Mass_tot1
a_3_z_E = (- N_motor * (2 + Operation_percentage[-1]) / 4 * np.cos(max(theta_reached_x, theta_reached_y) / np.sqrt(2)) * T2[7] + N_motor / 4 * np.cos(max(theta_reached_x, theta_reached_y) / np.sqrt(2)) * max(T_needed_x, T_needed_y) + Mass_tot1 * g) / Mass_tot1

V_x_E_end = V_x_E + a_3_x_E * t_stabilise_x
V_y_E_end = V_y_E + a_3_y_E * t_stabilise_y
V_Z_E_end = V_z_E + a_3_z_E * max(t_stabilise_x, t_stabilise_y)


print("Angle reached in the response time = ", theta_zero_vec_y * 57.3, theta_addition_y * 57.3)
print("Angular velocity to be counteracted = ", omega_tot_x, omega_tot_y)
print("Velocity of the craft = ", V_x_E, V_y_E)
print("Maximum positions = ", S_x_E) #y_max_E)
#print("Moments of inertia = ", I_xx, I_yy)

print("Thrusts = ", T_needed_x, T_needed_y)
print("Speed of craft relative to Earth = ", V_x_E_end, V_y_E_end, V_Z_E_end)
print("Accelerations relative to Earth", a_1_z_E, a_2_z_E, a_3_z_E)


"""    Gust load estimation    """








