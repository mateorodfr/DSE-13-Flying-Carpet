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
A_pipe = m_dot / rho_water / V_water
mat_perc = 1.05
D_pipe = np.sqrt(A_pipe / np.pi) * 2 * mat_perc

V_bat = concept.Mbat_concept / 1000 / n_HV_bat
print(D_pipe)
n_patch = 8

L_bat = 0.8
H_bat = 0.3
W_bat = V_bat / L_bat / H_bat
print(m_dot)
L_patch = L_bat / n_patch

S_cooling = 1.05 * ((n_patch + 1) * H_bat * W_bat + 2 * (L_bat * H_bat + L_bat * W_bat))
V_cooling = S_cooling * D_pipe
Perc_water = 0.2
M_water = V_cooling * Perc_water / mat_perc * n_HV_bat * rho_water
print(S_cooling * n_HV_bat)
print("The water mass in the pipes =", M_water, "kg")

S_mat_pipes = (D_pipe**2 - (D_pipe/mat_perc)**2) * np.pi
L_pipes = S_cooling * Perc_water / (D_pipe)
rho_pipes = 8730        # Brass applied
rho_plate = 2700
Plate_filling = 0.6
M_plate = V_cooling * (1 - Perc_water) * rho_plate * Plate_filling
M_pipes = S_mat_pipes * L_pipes * rho_pipes
M_tot = M_pipes + M_plate + M_water
print(M_tot)


