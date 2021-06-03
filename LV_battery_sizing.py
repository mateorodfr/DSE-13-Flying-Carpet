import numpy as np
import Parameters as pm

elect = pm.ElectronicsParameters(0)
concept = pm.ConceptParameters(0)

P_cons = elect.camera_power + elect.T_sens_power + elect.motor_controller_power + elect.VCU_power + elect.FC_power + elect.AMS_power + elect.SN_power
M_comp = elect.camera_mass * elect.camera_amount + elect.T_sens_mass * elect.T_sens_amount + elect.motor_controller_mass * elect.motor_controller_amount + elect.VCU_mass * (elect.VCU_amount - 1) + elect.FC_mass * (elect.FC_amount - 1) + elect.AMS_mass * (elect.AMS_amount - 1) + elect.SN_mass * elect.SN_amount

n_cyc = 20
t_asc_req = 80*n_cyc #time taken to ascend in s
t_des_req = 80*n_cyc #time taken to descend in s
t_hover_req = 445.5*n_cyc #time taken to hover in s
t_mission_req = (t_asc_req + t_hover_req + t_des_req) #total mission time

E_consumed = t_mission_req * P_cons / 3600

M_LV_bat = E_consumed / concept.battery.rhoE_battery / concept.battery.dod_battery / concept.battery.eff_battery / concept.battery.degradation
M_tot = M_LV_bat + M_comp
print(M_tot)
