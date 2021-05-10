import numpy as np
g=9.81
Maxpayload= 600
Missiont= 45/60

#Engine
Nengines= 4
BladeD=  np.arange(2.5,5.25,0.25)
RPM= 650
Torque= 1500
Power= 104000
CL= 1.1                                    #assumption NACA-2412
Nblade= 6
Sb= 0.1* BladeD**2 * Nblade
Vtip= 0.5*BladeD * 2* np.pi * RPM/60
Prop_eff= 0.85
M_eng= 49

#Reynolds number
rho= 1.225
V= Vtip / np.sqrt(2) #Change np.sqrt(2)
C= 0.1* BladeD
visc= 1.46*10**-5
Re= rho*V*C/visc

#Thrust per engine
T= Vtip**2 /6 * rho * CL * Sb *Prop_eff
print(T * Nengines /9.81)



#Battery stuff
Bat_E_dense=250 #Wh/kg
Bat_P_dense=500 #W/kg
Bat_eff=0.9
Eng_eff=0.85

#Power required
P_req= Power*Nengines *1.2
M_bat= P_req* Missiont/ Bat_E_dense /Bat_eff / Eng_eff
#print(M_bat)


Mass2= M_bat +Maxpayload + Nengines* M_eng + 1.5* Maxpayload

Bat_P= M_bat * Bat_P_dense

Bat_P_dense2= P_req/ M_bat



M_bat2= P_req / Bat_P_dense /Bat_eff /Eng_eff

M_total= 1.5*(Mass2+M_bat2)
print(M_bat2+M_bat)
print("Total mass", M_total)