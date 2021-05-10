import numpy as np
g=9.81
Maxpayload= 600

#Engine
Nengines= 4
BladeD= np.arange(2.5,5.25,0.25)
RPM= 650
Torque= 1500
Power= 204000
CL= 1.5                                    #assumption NACA-2412
Nblade= 5
Sb= 0.1* BladeD**2 * Nblade
Vtip= 0.5*BladeD * 2* np.pi * RPM/60

#Reynolds number
rho= 1.225
V= Vtip / np.sqrt(2) #Change np.sqrt(2)
C= 0.1* BladeD
visc= 1.46*10**-5

Re= rho*V*C/visc


T= Vtip**2 /6 * rho * CL * Sb
print(T * Nengines /9.81)




Thrust1= 600* 1.5  *9.81

P_T_ratio= 70000/1912.73529
P_W_ratio= 204000/49/9.81

Bat_E_dense=250 #Wh/kg
Bat_P_dense=750 #W/kg


