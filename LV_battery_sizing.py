import numpy as np
import Parameters as pm

elect = pm.ElectronicsParameters(0)
concept = pm.ConceptParameters(0)

V_cell = 3.7
n_cell = 4

P_cons = elect.pump_amount * elect.pump_power + elect.camera_power + elect.T_sens_power + elect.motor_controller_power + elect.VCU_power + elect.FC_power + elect.AMS_power + elect.SN_power
M_comp = elect.camera_mass * elect.camera_amount + elect.T_sens_mass * elect.T_sens_amount + elect.motor_controller_mass * elect.motor_controller_amount + elect.VCU_mass * (elect.VCU_amount - 1) + elect.FC_mass * (elect.FC_amount - 1) + elect.AMS_mass * (elect.AMS_amount - 1) + elect.SN_mass * elect.SN_amount

n_cyc = 20
t_asc_req = 80*n_cyc #time taken to ascend in s
t_des_req = 80*n_cyc #time taken to descend in s
t_hover_req = 445.5*n_cyc #time taken to hover in s
t_mission_req = (t_asc_req + t_hover_req + t_des_req) #total mission time

E_consumed = t_mission_req * P_cons / 3600      # Wh


V_LV = n_cell * V_cell

V_buck_12 = 12
V_buck_9 = 9
V_boost_24 = 24
f_s = 600 * 10**3
T_s = 1 / f_s
"""      The DCM (Discontinuous conduction mode) is selected for efficiency purposes and to meet the requirements
         the current ripple on the inductor has to be larger than twice the output power        """


P_out_12 = 2/3*(elect.VCU_power + elect.FC_power) + 1/2*elect.AMS_power + elect.pump_power * elect.pump_amount
I_out_12 = P_out_12/V_buck_12
L_buck_12_max = (np.sqrt((V_LV - V_buck_12) * T_s * 2 * P_out_12 / V_LV)/2/I_out_12)**2
L_buck_12 = np.arange(5*10**-9, L_buck_12_max, 5 * 10 ** -9)
D_12 = np.sqrt(2 * P_out_12 * L_buck_12 / (V_LV * (V_LV - V_buck_12) * T_s))
dI_L_12 = (V_LV - V_buck_12) * D_12 * T_s / L_buck_12       # The current ripple. Should be kept close to minimum to ensure
                                # that there is no overcurrent of the circuit, however, higher times of the zero current
                                # are benificial for efficiency

d_V_out = 0.05

C_12 = 1 / (2 * d_V_out) * (dI_L_12 - I_out_12) * (L_buck_12 * (dI_L_12 - I_out_12) / V_buck_12 + D_12 * T_s - L_buck_12 * I_out_12 / (V_LV - V_buck_12))
d12 = D_12 * (V_LV/V_buck_12 - 1)
t_zero_12 = T_s - d12 * T_s - D_12 * T_s
I_rms_L_12 = np.sqrt((D_12 + d12) * dI_L_12**2 / 3)
I_rms_sw_12 = np.sqrt(D_12 * dI_L_12**2 / 3)
I_rms_d_12 = np.sqrt(d12 * dI_L_12**2 / 3)
I_avg_d_12 = d12 * dI_L_12 / 2
I_rms_C_12 = np.sqrt((D_12 + d12) * (dI_L_12**2/3 - dI_L_12 * I_out_12) + I_out_12**2)
R_L_12 = 0.005
R_sw_12 = 0.0047
R_d_12 = 0.003
V_d_on_12 = 0.7
R_C_12 = 0.025
P_loss_12 = I_avg_d_12 * V_d_on_12 + I_rms_d_12**2 * R_d_12 + R_sw_12 * I_rms_sw_12**2 + R_L_12 * I_rms_L_12**2 + R_C_12 * I_rms_C_12**2
eff_buck_12 = (P_out_12 - P_loss_12) / P_out_12


print(eff_buck_12)

"""    Rest of the VCU and AMS systems work on 24V, thus Boost converter is applied     """

d_V_out = 0.05

P_out_boost_24 = 1/3*elect.VCU_power + 1/2*elect.AMS_power

I_in_24 = P_out_boost_24 / V_LV
I_out_24 = P_out_boost_24 / V_boost_24
L_boost_24_max = (2 * P_out_boost_24 * T_s * (V_boost_24 - V_LV) / V_boost_24) / (4 * I_in_24**2)
L_boost_24 = np.arange(20 * 10 ** -9, L_boost_24_max, 40 * 10 ** -9)
dI_L_24 = np.sqrt(2 * P_out_boost_24 * T_s * (V_boost_24 - V_LV) / (V_boost_24 * L_boost_24))
D_24 = dI_L_24 * L_boost_24 / V_LV / T_s
d24 = 2 * I_out_24 / dI_L_24

C_24 = 1 / (2 * d_V_out) * L_boost_24 * (dI_L_24 * I_out_24)**2/(V_boost_24 - V_LV)
I_rms_L_24 = np.sqrt((D_24+d24) * dI_L_24**2 / 3)
I_rms_sw_24 = np.sqrt(D_24 * dI_L_24**2 / 3)
I_rms_d_24 = np.sqrt(d24 * dI_L_24**2 / 3)
I_avg_d_24 = d24 * dI_L_24 / 2
I_rms_C_24 = np.sqrt((d24/3 - d24**2/4) * dI_L_24**2)
R_L_24 = 0.005
R_sw_24 = 0.0047
R_d_24 = 0.003
V_d_on_24 = 0.7
R_C_24 = 0.025
P_loss_24 = I_avg_d_24 * V_d_on_24 + I_rms_d_24**2 * R_d_24 + R_sw_24 * I_rms_sw_24**2 + R_L_24 * I_rms_L_24**2 + R_C_24 * I_rms_C_24**2
eff_boost_24 = (P_out_boost_24 - P_loss_24) / P_out_boost_24

print(eff_boost_24)


P_out_9 = elect.SN_power * elect.SN_amount + 1/3 * elect.FC_power + elect.camera_power * elect.camera_power
I_out_9 = P_out_9/V_buck_9
L_buck_9_max = (np.sqrt((V_LV - V_buck_9) * T_s * 2 * P_out_9 / V_LV)/2/I_out_9)**2
L_buck_9 = np.arange(5*10**-9, L_buck_9_max, 5 * 10 ** -9)
D_9 = np.sqrt(2 * P_out_9 * L_buck_9 / (V_LV * (V_LV - V_buck_9) * T_s))
dI_L_9 = (V_LV - V_buck_9) * D_9 * T_s / L_buck_9       # The current ripple. Should be kept close to minimum to ensure
                                # that there is no overcurrent of the circuit, however, higher times of the zero current
                                # are benificial for efficiency

C_9 = 1 / (2 * d_V_out) * (dI_L_9 - I_out_9) * (L_buck_9 * (dI_L_9 - I_out_9) / V_buck_9 + D_9 * T_s - L_buck_9 * I_out_9 / (V_LV - V_buck_9))
d9 = D_9 * (V_LV/V_buck_9 - 1)
t_zero_9 = T_s - d9 * T_s - D_9 * T_s
I_rms_L_9 = np.sqrt((D_9 + d9) * dI_L_9**2 / 3)
I_rms_sw_9 = np.sqrt(D_9 * dI_L_9**2 / 3)
I_rms_d_9 = np.sqrt(d9 * dI_L_9**2 / 3)
I_avg_d_9 = d9 * dI_L_9 / 2
I_rms_C_9 = np.sqrt((D_9 + d9) * (dI_L_9**2/3 - dI_L_9 * I_out_9) + I_out_9**2)
R_L_9 = 0.005
R_sw_9 = 0.0047
R_d_9 = 0.003
V_d_on_9 = 0.7
R_C_9 = 0.025
P_loss_9 = I_avg_d_9 * V_d_on_9 + I_rms_d_9**2 * R_d_9 + R_sw_9 * I_rms_sw_9**2 + R_L_9 * I_rms_L_9**2 + R_C_9 * I_rms_C_9**2
eff_buck_9 = (P_out_9 - P_loss_9) / P_out_9

print(eff_buck_9)

M_LV_bat = E_consumed / concept.battery.rhoE_battery / concept.battery.dod_battery / concept.battery.eff_battery  / ((eff_buck_9[-1] * P_out_9 + eff_buck_12[-1] * P_out_12 + eff_boost_24[-1] * P_out_boost_24)/P_cons) / elect.PF_electronics
M_tot = M_LV_bat + M_comp
print(P_cons/M_LV_bat)

print(M_tot)

