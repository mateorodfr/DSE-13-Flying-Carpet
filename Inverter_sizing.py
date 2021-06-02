import numpy as np
import Parameters as pm

concept = pm.ConceptParameters(0)

n_cyc = 1
t_asc_req = 80*n_cyc #time taken to ascend in s
t_des_req = 80*n_cyc #time taken to descend in s
t_hover_req = 445.5*n_cyc #time taken to hover in s
t_mission_req = (t_asc_req + t_hover_req + t_des_req) #total mission time

V_cell = 3.7
V_engine = np.arange(450, 851, 10)
P_engine = np.arange(0, 102000, 2500)
n_cell = 250
Cap_dens = concept.battery.rhoE_battery / V_cell
M_bat_it = 471     # Got from mass estimation repeated
P_cont_tot = 434210
#P_surge = 1.5 * P_cont_tot
m_row = M_bat_it/n_cell
Cap_tot = m_row * Cap_dens
#I_bat = Cap_tot / t_mission_req * 3600
Inv_size = P_cont_tot * concept.battery.degradation * concept.battery.eff_battery / concept.motor.N_motor# VA
Bat_current = P_cont_tot / (V_cell * n_cell)

"""     Both the delta and star configurations have to be applied to the motor in order to reduce the inrush/surge 
        current of the motor and to have a safety mechanism/interlock to avoid burning of batteries. Star configuration
        is used first as the voltage and the current are both lower and then the delta configuration takes over to
        increase the power output for the motor     """

R_engine_wire = V_engine[-1]**2 / P_engine[-1]
print(R_engine_wire)
V_line_inv = np.sqrt(Inv_size * R_engine_wire)
print(V_line_inv)
l_line_inv = 2.5
I_line_inv = Inv_size / V_line_inv
rho_wire = 1.724 * 10 ** -8

#P_loss_inv = I_line_inv ** 2 * R_wire_inv
#P_max = np.max(P_loss_inv)

S_wire = (R_engine_wire / rho_wire / l_line_inv)**(-1)
d_wire = np.sqrt(S_wire / np.pi) * 2
V_wire_inv = l_line_inv * S_wire
print(d_wire)
rho_copper = 8940
m_wire = rho_copper * V_wire_inv

I_delta_coil = I_line_inv / np.sqrt(3)
R_delta_coil = V_line_inv / I_delta_coil

I_line_inv_startup = V_engine[0] / R_engine_wire
V_star_coil = V_engine[0] / np.sqrt(3)
R_star = V_star_coil / I_line_inv_startup

Q_wire_max = Inv_size * t_mission_req
c = 385
dT = Q_wire_max / c / np.max(m_wire)

print(I_line_inv)



