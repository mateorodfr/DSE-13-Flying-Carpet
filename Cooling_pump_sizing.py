import Parameters as pm

concept = pm.ConceptParameters(0)

C_Li = 3582     # J/kg K
C_Epoxy = 1110  # J/kg K
C_V_O = 490     # J/kg K
n_Li = 1/6
n_V_0 = 1/6
n_Epoxy = 2/3
n_HV_bat = 8

t_asc_req = 60
t_des_req = 60
t_hover_req = 445.5
t_mission = (t_asc_req + t_hover_req + t_des_req)

C_tot = n_V_0 * C_V_O + n_Li * C_Li + n_Epoxy * C_Epoxy
P_loss = concept.Mbat_concept * concept.battery.rhoP_battery * concept.battery.eff_battery
Q_loss = P_loss * t_mission
dT_per_HV = Q_loss / (C_tot * concept.Mbat_concept / n_HV_bat)
print(P_loss)

