import numpy as np
#waddup
g = 9.81
Maxpayload = 600
Missiont = 45/60


#this is my new comment
#I dont want to push this comment though
#Engine
Nengines = 8
R_engine = 0.418 / 2
BladeD = np.arange(2.5, 3, 0.1)
#RPM = 1050
#omega = RPM /60 * 2 * np.pi
#freq = omega / 2 / np.pi
Torque = 1000  # TBD
#Power = 104000    # W
CL = 0.9   #assumption NACA-2412
Nblade = 6
Sb = 0.1* BladeD**2 * Nblade
t_blade = 0.01 * BladeD / 2

Prop_eff= 0.85
Eng_eff=0.85
M_eng= 49
rho_blade = 660 #kg/m^3
M_blades = Sb * t_blade * rho_blade
I_tot = 1/2 * M_eng * R_engine**2 + 1/3 * M_blades * (BladeD*BladeD)/4
omega = Torque/I_tot

Vtip= 0.5*BladeD * 2 * omega
P_eng = Torque * omega
print("P_eng = ", P_eng)

#Reynolds number
rho= 1.225
V= Vtip / np.sqrt(2)    #Change np.sqrt(2)
C= 0.1* BladeD
visc= 1.46*10**-5
Re= rho*V*C/visc

#Thrust per engine
T = Vtip**2 /6 * rho * CL * Sb *Prop_eff
print("Mass to be carried by engines = ", T * Nengines /9.81)

#Battery stuff
Bat_E_dense=250 #Wh/kg
#Bat_P_dense=500 #W/kg
Bat_eff=0.9


#Power required
P_req = P_eng * Nengines / Eng_eff / Prop_eff
M_bat = P_req * Missiont / Bat_E_dense / Bat_eff / Eng_eff
#print(M_bat)

Mass_tot = 1.2 * (M_bat + Maxpayload + Nengines * M_eng + M_blades)


print("Total mass = ", Mass_tot)
print("Safety margin = ", T/Mass_tot)