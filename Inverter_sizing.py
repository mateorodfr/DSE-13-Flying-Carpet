import numpy as np
import Parameters as pm

concept = pm.ConceptParameters(0)

n_cyc = 1
t_asc_req = 80*n_cyc #time taken to ascend in s
t_des_req = 80*n_cyc #time taken to descend in s
t_hover_req = 445.5*n_cyc #time taken to hover in s
t_mission_req = (t_asc_req + t_hover_req + t_des_req) #total mission time

V_cell = 3.7
n_cell = 200
Cap_dens = concept.battery.rhoE_battery / V_cell
M_bat_it = 471     # Got from mass estimation repeated
P_cont_tot = 434210
#P_surge = 1.5 * P_cont_tot
m_row = M_bat_it/n_cell
Cap_tot = m_row * Cap_dens
#I_bat = Cap_tot / t_mission_req * 3600
Inv_size = P_cont_tot * concept.battery.degradation * concept.battery.eff_battery # VA
Bat_current = P_cont_tot / (V_cell * n_cell)

l_line = 4
L_tot = l_line * concept.motor.N_motor
d_wire = 0.006
S_wire = (d_wire/2)**2 * np.pi
rho_wire = 1.68 * 10 ** -8
R_wire = rho_wire * L_tot / S_wire
P_loss = R_wire * Bat_current ** 2
V_wire = S_wire * L_tot
rho_copper = 8940
m_wire = rho_copper * V_wire

Q_wire = P_loss * t_mission_req
c = 385
dT = Q_wire / c / m_wire

print(Bat_current, P_loss, dT)



