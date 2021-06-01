import numpy as np
import Parameters as pm
import matplotlib.pyplot as plt
from math import *


def motion(r0,v0,a0,t):

    r1 = r0 + v0*t + 0.5*a0*t*t
    v1 = a0*t

    return r1,v1



concept = pm.ConceptParameters(0)


#Operational times (arbitrary)
n_cyc = 1
t_asc_req = 80*n_cyc #time taken to ascend in s
t_des_req = 80*n_cyc #time taken to descend in s
t_hover_req = 445.5*n_cyc #time taken to hover in s
t_mission_req = (t_asc_req + t_hover_req + t_des_req) #total mission time

#
dt = 0.01
h_max = 400


#Propulsive times

#OEI acc needed
a_OEI = 1.6 * concept.physics.g
a_roll = 0.4 * concept.physics.g

#Takeoff initial state
a0 = 0.1*concept.physics.g #Takeoff acceleration
r0 = 0  #Initial position
v0 = 0  #Initial velocity

#Descent initial state
a2 = -0.1*concept.physics.g
r2 = h_max
v2 = 0

def motionTime(r0,v0,a0,r_max,t_req,dt):
    t_tot = 10000
    t0 = 0
    while np.abs(t_tot) > t_req:
    #Manoevure time
    #Ascent at constant speed
        t0 += dt
        r1,v1 = motion(r0,v0,a0,t0)
        t1 = np.abs(r_max-2*r1)/v1
        t_tot = 2*t0 + t1
    return t_tot,t0

t_asc,t_asc_acc = motionTime(r0,v0,a0,h_max,t_asc_req,dt)
t_asc_noa = t_asc-2*t_asc_acc
t_hover = t_hover_req
t_des = t_asc
t_des_acc = t_asc_acc
t_des_noa = t_asc_noa

t_mission = t_asc + t_hover + t_des
t = np.arange(0,t_mission+dt,dt)

# Assumptions for now

M_tot0 = 1200
M_tot1 = concept.Mtot_concept
CL = concept.propeller.CL_prop
Sb = concept.propeller.S_prop
Dblade = concept.propeller.D_prop
torque = concept.motor.Torque

iterations=0
while abs(1-M_tot1/M_tot0) > 0.01:
    progress = 0
    state = []
    T_req_eng=[]
    for i in range(len(t)):


        if np.round((i/len(t))*100) > progress:
            progress = np.round((i/len(t))*100)
            #print(f'Progress: {progress} %')
        if t[i] < t_asc:
            if t[i] < t_asc_acc:
                ai = 0.1*concept.physics.g
                vi = ai*t[i]
                ri = 0.5*ai*t[i]*t[i]
                ti = t[i]
                state.append([ai,vi,ri,ti])
                T_req_eng.append((M_tot1-500)*(ai + concept.physics.g)/concept.motor.N_motor /concept.propeller.eff_prop)

            elif t_asc_acc <= t[i] < t_asc_noa+t_asc_acc:
                aj = 0
                vj = aj*(t[i]-t_asc_acc)+vi
                rj = ri + vi*(t[i]-t_asc_acc) + 0
                tj = t[i]
                state.append([aj,vj,rj,tj])
                T_req_eng.append((M_tot1 - 500) * (
                            aj + concept.physics.g) / concept.motor.N_motor / concept.propeller.eff_prop)

            elif t_asc_noa <= t[i] < t_asc:
                ak = -0.1*concept.physics.g
                vk = vj + ak*(t[i]-t_asc_noa-t_asc_acc)
                rk = rj + vj*(t[i]-t_asc_noa-t_asc_acc) + 0.5*ak*(t[i]-t_asc_noa-t_asc_acc)**2
                tk = t[i]
                state.append([ak,vk,rk,tk])
                T_req_eng.append((M_tot1 - 500) * (
                            ak + concept.physics.g) / concept.motor.N_motor / concept.propeller.eff_prop)

        elif t_asc <= t[i] < t_hover+t_asc:
            al = 0
            vl = 0
            rl = rk
            tl = t[i]
            state.append([al,vl,rl,tl])
            T_req_eng.append((M_tot1) * (
                        al + concept.physics.g) / concept.motor.N_motor / concept.propeller.eff_prop)

        elif t_hover+t_asc <= t[i] < t_mission:
            if t[i] < t_hover+t_asc+t_des_acc:
                ac = -0.1*concept.physics.g
                vc = vl + ac*(t[i]-t_hover-t_asc)
                rc = rl + 0.5*ac*(t[i]-t_hover-t_asc)**2
                tc = t[i]
                state.append([ac,vc,rc,tc])
                T_req_eng.append((M_tot1) * (
                        ac + concept.physics.g) / concept.motor.N_motor / concept.propeller.eff_prop)

            elif t_hover+t_asc+t_des_acc <= t[i] < t_hover+t_asc+t_des_acc+t_des_noa:
                av = 0
                vv = av*(t[i]-t_hover-t_asc-t_des_acc)+vc
                rv = rc + vc*(t[i]-t_hover-t_asc-t_des_acc) + 0
                tv = t[i]
                state.append([av,vv,rv,tv])
                T_req_eng.append((M_tot1) * (
                        av + concept.physics.g) / concept.motor.N_motor / concept.propeller.eff_prop)

            elif t_hover+t_asc+t_des_acc+t_des_noa <= t[i] < t_mission:
                ab = 0.1*concept.physics.g
                vb = vv + ab*(t[i]-t_hover-t_asc-t_des_acc-t_des_noa)
                rb = rv + vv*(t[i]-t_hover-t_asc-t_des_acc-t_des_noa) + 0.5*ab*(t[i]-t_hover-t_asc-t_des_acc-t_des_noa)**2
                tb = t[i]
                state.append([ab,vb,rb,tb])
                T_req_eng.append((M_tot1) * (
                        ab + concept.physics.g) / concept.motor.N_motor / concept.propeller.eff_prop)

    # r0 = v0 = a0 = np.zeros(len(t))
    # a00 = a0[np.where((t>=0)&(t<t_asc_acc))]+0.1*concept.physics.g
    # a01 = a0[np.where((t>=t_asc_acc)&(t<t_asc_acc+t_asc_noa))]
    # a02 = a0[np.where((t>=t_asc_acc+t_asc_noa)&(t<t_asc))]-0.1*concept.physics.g
    # a10 = a0[np.where((t>=t_asc)&(t<t_asc+t_hover))]
    # a20 = a0[np.where((t>=t_asc+t_hover)&(t<t_asc+t_hover+t_des_acc))]-0.1*concept.physics.g
    # a21 = a0[np.where((t>=t_asc+t_hover+t_des_acc)&(t<t_asc+t_hover+t_des_acc+t_des_noa))]
    # a22 = a0[np.where((t>=t_asc+t_hover+t_des_acc+t_des_noa)&(t<t_mission))]+0.1*concept.physics.g
    # a_mission = np.append(np.append(np.append(np.append(np.append(np.append(a00,a01),a02),a10),a20),a21),a22)


    # print(ai,vi,ri,ti)
    # print(aj,vj,rj,tj)
    # print(ak,vk,rk,tk)
    # print(al,vl,rl,tl)
    # print(ac,vc,rc,tc)
    # print(av,vv,rv,tv)
    # print(ab,vb,rb,tb)

    state = np.array(state)





    #T_req_eng = M_tot0*(state[:,0]+concept.physics.g)/concept.motor.N_motor /concept.propeller.eff_prop
    T_req_eng= np.array(T_req_eng)
    #P_req_eng= np.sqrt(T_req_eng/concept.physics.rho0/CL/Sb -0.25* state[:,1]**2) *3* torque/ np.pi**2/Dblade /concept.motor.eff_motor

    P_req_eng = np.sqrt((2*T_req_eng**3)/(concept.physics.rho0*np.pi*(concept.propeller.D_prop/concept.propeller.eff_prop)**2))/(concept.motor.eff_motor) * (1+ concept.battery.loss_factor) / concept.battery.eff_inverter / concept.battery.degradation
    P_req_eng_OEI = sqrt((2 * (M_tot1*a_OEI/2 + M_tot1*a_roll/4) ** 3) / (concept.physics.rho0 * np.pi * (concept.propeller.D_prop / concept.propeller.eff_prop) ** 2)) / (concept.motor.eff_motor) / concept.battery.loss_factor / concept.battery.eff_inverter

    E_req = np.sum((concept.motor.N_motor*P_req_eng)*dt) / concept.battery.eff_battery / concept.battery.dod_battery * (1 + concept.battery.loss_factor) / concept.battery.eff_inverter / concept.battery.degradation
    Mbat = E_req/ (concept.battery.rhoE_battery *3600)


    if Mbat < np.max(P_req_eng) * concept.motor.N_motor / concept.battery.rhoP_battery:
        Mbat = np.max(P_req_eng) * concept.motor.N_motor / concept.battery.rhoP_battery

    # if Mbat < P_req_eng_OEI / concept.battery.rhoP_battery:
    #     Mbat = P_req_eng_OEI / concept.battery.rhoP_battery



    M_tot0= M_tot1


    M_motor= concept.motor.M_motor * concept.motor.N_motor
    M_propeller= concept.propeller.M_blades* concept.motor.N_motor
    M_payload= concept.Mpay_concept
    M_struct_coeff= 1.2

    M_tot1= (Mbat + M_motor + M_propeller +M_payload)*M_struct_coeff
    iterations +=1

print("Final mass",M_tot1, "After", iterations, "Iterations")
print("Battery mass",Mbat)
print("---------------")
print("Required power", np.max(P_req_eng)*concept.motor.N_motor )
print("OEI power", np.max(P_req_eng_OEI))
print("Required power density",np.max(P_req_eng*concept.motor.N_motor)/Mbat, P_req_eng_OEI/Mbat)
print("Power density battery", concept.battery.rhoP_battery)
print("---------------")
print("Required energy", E_req)
print("Required energy density", E_req/Mbat)
print("Energy density battery", concept.battery.rhoE_battery*3600)