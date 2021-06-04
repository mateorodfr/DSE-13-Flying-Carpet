import numpy as np
import Parameters as pm

concept = pm.ConceptParameters(0)

C_Li = 3582     # J/kg K
C_Epoxy = 1110  # J/kg K
C_V_O = 490     # J/kg K
n_Li = 1/6
n_V_0 = 1/6
n_Epoxy = 2/3
n_HV_bat = 4

rho_pipes = 8730        # Copper-Nickel 70/30 applied
rho_plate = 2700        # Al applied
Plate_filling = 1
Perc_water = 0.2

t_asc_req = 60
t_des_req = 60
t_hover_req = 445.5
t_mission = (t_asc_req + t_hover_req + t_des_req)

C_tot = n_V_0 * C_V_O + n_Li * C_Li + n_Epoxy * C_Epoxy
P_loss = 1.05 * concept.Mbat_concept * concept.battery.rhoP_battery * concept.battery.eff_battery / n_HV_bat
Q_loss = P_loss * t_mission
dT_per_HV = Q_loss / (C_tot * concept.Mbat_concept)
#print(dT_per_HV)

eff_cooling = 0.7       # Due the brass/copper presence and the water friction with the pipes
T_min = 20  # Celsius
T_max = 50  # Celsius
Cp = 4186   # J/kg K
dT = T_max - T_min
m_dot = P_loss/(dT*Cp) / eff_cooling     # kg/S
rho_water = 1000    # kg/m^3
V_water = 4.6   # m/s
mat_perc = 1.05 ** 2
A_pipe = m_dot / rho_water / V_water * mat_perc
D_pipe = np.sqrt(A_pipe / np.pi) * 2

V_bat = concept.Mbat_concept / 1000 / n_HV_bat
print(D_pipe)
n_patch = 8

L_bat = 0.8
H_bat = 0.3
W_bat = V_bat / L_bat / H_bat
print(m_dot)
L_patch = L_bat / n_patch

h = 350
A_radiator = P_loss / h / dT
print(A_radiator)

t_cool = 0.00035
gap = 2 * t_cool
L_cool = 0.6          # Length of the small plate
W_cool = 0.08       # Depth of the small plate
H_cool = A_radiator / (W_cool * L_cool) * (t_cool + gap) * 1.1   # Height of the radiator

V_rad = L_cool * t_cool * A_radiator / (W_cool * L_cool) * W_cool
extra_meters = 3

D_cool_pipes = W_cool/1.2
A_cool_pipes = np.pi * (((D_cool_pipes * mat_perc)**2 - (D_cool_pipes)**2)/4)
L_cool_pipes = H_cool / 2 / D_cool_pipes * L_cool
V_cool_pipes = L_cool_pipes * A_cool_pipes
m_cool_pipes = V_cool_pipes * rho_plate

A_cooling = 1.05 * ((n_patch + 1) * H_bat * W_bat + 2 * (L_bat * H_bat + L_bat * W_bat))
V_cooling = A_cooling * D_pipe
M_water = (V_cooling * Perc_water / mat_perc + A_pipe / mat_perc * extra_meters + V_cool_pipes / mat_perc) * n_HV_bat * rho_water
print(A_cooling * n_HV_bat)
print("The water mass in the pipes =", M_water, "kg")

A_mat_pipes = A_pipe - A_pipe/mat_perc
L_pipes = (V_cooling * Perc_water * (1 - 1/mat_perc))/A_mat_pipes + extra_meters
M_plate = V_cooling * (1 - Perc_water) * rho_plate * Plate_filling * n_HV_bat
M_pipes = A_mat_pipes * L_pipes * rho_pipes * n_HV_bat
M_tot = M_pipes + M_plate + M_water
print(M_plate)

m_rad = (V_rad * rho_plate + m_cool_pipes) * 1.1
M_pump = 1.3
M_final = m_rad + M_tot + n_HV_bat * M_pump
print(M_final)

