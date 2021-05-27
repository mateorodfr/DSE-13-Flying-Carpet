import Parameters as pm
import numpy as np

#Import parameters
concept = pm.ConceptParameters(0) #Concept
physics = pm.PhysicalParameters() #Physical constants

P_eng = np.arange(0,concept.motor.P_max+1,concept.motor.P_max/1000).astype(np.float32) #Single motor power range
T_eng = ((np.pi/2)*concept.propeller.D_prop**2*physics.rho0*(P_eng*concept.motor.eff_motor*concept.propeller.eff_prop)**2)**(1/3) #Single motor thrust range


#Concept Masses from initial sizing
m_tot = concept.Mtot_concept
m_bat = concept.Mbat_concept
m_pay = concept.Mpay_concept

#Thurst over weight ration range
TW = (concept.motor.N_motor*T_eng)/(m_tot*physics.g)

#Indices of power setting for operations
idx_ops = [np.abs(TW - 1.1).argmin(),np.abs(TW - 1.0).argmin(),np.abs(TW - 0.9).argmin()]
P_ops = np.array([P_eng[i] for i in idx_ops])
T_ops = np.array([T_eng[i] for i in idx_ops])
TW_ops = np.array([TW[i] for i in idx_ops])
Hover_efficiency = (T_ops*1000)/(P_ops*physics.g)

E_tot = m_bat*concept.battery.rhoE_battery*3600*concept.battery.eff_battery #Energy stored in joules
t_ops = E_tot/(concept.motor.N_motor*P_ops)

isPrint = True
if isPrint:
    print('\nOperation time [min]: ', t_ops/60 ,'\nOperational Power [kW]: ' , P_ops/1000 , '\nOperational Thrust-over-Weight ratio [-]:' , TW_ops, '\nHover efficiency [kg/W]: ', Hover_efficiency)

