import numpy as np
import Parameters as pm
import matplotlib.pyplot as plt



def motion(r0,v0,a0,t):

    r1 = r0 + v0*t + 0.5*a0*t*t
    v1 = a0*t

    return r1,v1



concept = pm.ConceptParameters(0)


#Operational times (arbitrary)
n_cyc = 1
t_asc_req = 60*n_cyc #time taken to ascend in s
t_des_req = 60*n_cyc #time taken to descend in s
t_hover_req = 11*60*n_cyc #time taken to hover in s
t_mission_req = (t_asc_req + t_hover_req + t_des_req) #total mission time

#
dt = 0.01
t = np.arange(0,t_mission_req+dt,dt)
h_max = 400


#Propulsive times

#Takeoff initial state
a0 = 0.1*concept.physics.g #Takeoff acceleration
r0 = 0  #Initial position
v0 = 0  #Initial velocity

#Descent initial state
a2 = -0.1*concept.physics.g
r2 = h_max
v2 = 0

t_tot = 1000
t0 = 0



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



progress = 0
state = []
for i in range(len(t)):


    if np.round((i/len(t))*100) > progress:
        progress = np.round((i/len(t))*100)
        print(f'Progress: {progress} %')
    if t[i] < t_asc:
        if t[i] < t_asc_acc:
            ai = 0.1*concept.physics.g
            vi = ai*t[i]
            ri = 0.5*ai*t[i]*t[i]
            ti = t[i]
            state.append([ai,vi,ri,ti])

        elif t_asc_acc <= t[i] < t_asc_noa+t_asc_acc:
            aj = 0
            vj = aj*(t[i]-t_asc_acc)+vi
            rj = ri + vi*(t[i]-t_asc_acc) + 0
            tj = t[i]
            state.append([aj,vj,rj,tj])

        elif t_asc_noa <= t[i] < t_asc:
            ak = -0.1*concept.physics.g
            vk = vj + ak*(t[i]-t_asc_noa-t_asc_acc)
            rk = rj + vj*(t[i]-t_asc_noa-t_asc_acc) + 0.5*ak*(t[i]-t_asc_noa-t_asc_acc)**2
            tk = t[i]
            state.append([ak,vk,rk,tk])
    elif t_asc <= t[i] < t_hover+t_asc:
        al = 0
        vl = 0
        rl = rk
        tl = t[i]
        state.append([al,vl,rl,tl])

    elif t_hover+t_asc <= t[i] < t_mission:
        if t[i] < t_hover+t_asc+t_des_acc:
            ac = -0.1*concept.physics.g
            vc = vl +  ac*(t[i]-t_hover-t_asc)
            rc = rl + 0.5*ac*(t[i]-t_hover-t_asc)**2
            tc = t[i]
            state.append([ac,vc,rc,tc])

        elif t_hover+t_asc+t_des_acc <= t[i] < t_hover+t_asc+t_des_acc+t_des_noa:
            av = 0
            vv = av*(t[i]-t_hover-t_asc-t_des_acc)+vc
            rv = rc + vc*(t[i]-t_hover-t_asc-t_des_acc) + 0
            tv = t[i]
            state.append([av,vv,rv,tv])

        elif t_hover+t_asc+t_des_acc+t_des_noa <= t[i] < t_mission:
            ab = 0.1*concept.physics.g
            vb = vv + ab*(t[i]-t_hover-t_asc-t_des_acc-t_des_noa)
            rb = rv + vv*(t[i]-t_hover-t_asc-t_des_acc-t_des_noa) + 0.5*ab*(t[i]-t_hover-t_asc-t_des_acc-t_des_noa)**2
            tb = t[i]
            state.append([ab,vb,rb,tb])

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
T_req_eng = concept.Mtot_concept*(state[:,0]+concept.physics.g)/concept.motor.N_motor
P_req_eng = np.sqrt((2*T_req_eng**3)/(concept.physics.rho0*np.pi*concept.propeller.D_prop**2))/(concept.propeller.eff_prop*concept.motor.eff_motor)


E = np.sum((concept.motor.N_motor*P_req_eng)*dt)/concept.battery.eff_battery/concept.battery.dod_battery
mbat = E/(concept.battery.rhoE_battery*3600)
print(np.max(P_req_eng*concept.motor.N_motor)/concept.Mbat_concept)
print(mbat,concept.Mbat_concept/0.7)





plt.subplot(3,2,1)
plt.plot(state[:,-1],state[:,0])
plt.title('Acceleration over time for one cycle')
plt.xlabel('Time [s]')
plt.ylabel(r'Acceleration [ms$^{-2}$]')
plt.subplot(3,2,2)
plt.plot(state[:,-1],state[:,1])
plt.title('Vertical velocity over time for one cycle')
plt.xlabel('Time [s]')
plt.ylabel(r'Velocity [ms$^{-1}$]')
plt.subplot(3,2,3)
plt.plot(state[:,-1],state[:,2])
plt.title('Altitude over time for one cycle')
plt.xlabel('Time [s]')
plt.ylabel(r'Altitude [m]')
plt.subplot(3,2,4)
plt.plot(state[:,-1],T_req_eng)
plt.title('Thrust over time for one cycle per engine')
plt.xlabel('Time [s]')
plt.ylabel(r'Thrust [N]')
plt.subplot(3,2,5)
plt.plot(state[:,-1],P_req_eng)
plt.title('Power over time for one cycle per engine')
plt.xlabel('Time [s]')
plt.ylabel(r'Power [W]')
plt.show()




a1 = 0

