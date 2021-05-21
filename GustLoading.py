import numpy as np
import Parameters as pm
import matplotlib.pyplot as plt
import MomentInertia as mi


"""Define Module Specific parameters"""

#Gust load constants
Ui = 5 #wind speed at sensor height in [m/s]
z0 = 3 #surface roughness constant -> 3m for urban centers [m]
zi = 10 #Anemometer height when measuring wind spped of ground [m]
zf = 400 #Maximum operational height [m]

dz = 0.5 #height step [m]
z = np.arange(zi,zf+dz,dz) #operational height range [m]
Uf = Ui*(np.log(z/z0)/np.log(zi/z0)) #Gust speed at given height [m/s]
Umax = Uf[-1] #Max gust load

""" Import General Parameters"""

#Import subsystem parameters
cabin = pm.CabinParameters(0)
motor = pm.MotorParameters(0)
propeller = pm.PropellerParameters(0)
concept = pm.ConceptParameters(0)

#System parameters
#all these variable sneed to be manually changed to perform a case specific analysis
Cd = 1.2 #CD of squre at reynolds: 6e5
rho = 1.255 #[kg/m3] density at sea level
H_engine = 0.3
m = 1455 #curent final iterated mass
Mpayload = 600 #Current max payload mass
mOEW = m - Mpayload - motor.N_motor*(motor.M_motor+propeller.M_blades)


#Payload/Person dimensions
L_person = cabin.L_person
W_person = cabin.W_person
H_person = cabin.H_person

#Cabin dimensions
L_cabin = cabin.L_cabin #Cabin length
W_cabin = cabin.W_cabin #Cabin width
H_cabin = cabin.H_cabin #Cabin height
Dim_cabin = np.array([L_cabin,W_cabin,H_cabin]) #Array to store all cabin dimensions

#Cabin surface Areas
Syz = cabin.S_cabin[2] #yz plane
Sxz = cabin.S_cabin[1] #xz plane
Sxy = cabin.S_cabin[0] #xy plane
S = np.array([Syz,Sxz,Sxy]) #Array to store surface areas in x,y,z direction


""" Plot gusts"""
plot = False
if plot:
    plt.plot(z,Uf)
    plt.title('Gust speed at varying altitudes')
    plt.xlabel('Altitude: h [m]')
    plt.ylabel(r'Wind speed: U [ms$^{-1}$]')
    plt.savefig("figures/gustload")

"""Compute Moment of Inertias"""
I = mi.MMOI(motor,propeller,cabin,concept)

"""Compute characteristics"""

#Maximum gust force
Dgust = 0.5*rho*Umax**2*S*Cd

#Translational accelerations gust point load applied at cg
a = Dgust/m

#Rotational accelerations due to distributed gust load.
"""
Assuming the same magnitude of gust but distributed in a triangle.
This will cause a translational acceleration equal to ax,ay,az calculated above.
The distributed load will cause a torque.
The load is simplified as a point load applied at 1/3 of the load.
Assuming the load goes from 0 to w0 over the cabin
This means the load is applied at 1/3 Lcabin, Wcabin, Hcabin
"""
d = (1/6)*Dim_cabin #Moment arm
M  = Dgust*d #Torque
alpha = M/I #Angular accelerations



"""Print Characteristics"""
print(
    f'The gust load analysis yields: '
    f'\n\tThe moment of inertias in x,y,z'
    f'\n\t\t Ixx: {np.round(I[0])} [kgm2]'
    f'\n\t\t Iyy: {np.round(I[1])} [kgm2]'
    f'\n\t\t Izz: {np.round(I[2])} [kgm2]'
    f'\n\tThe maximum gust speed in x,y,z:'
    f'\n\t\t Umax: {np.round(Umax,2)} [m/s]'
    f'\n\tThe maximum distrubance force due to gust load in x,y,z:'
    f'\n\t\t Fx, Fy, Fz: {np.round(Dgust,2)} [N]'
    f'\n\tThe maximum translational accelerations due to disturbance in x,y,z'
    f'\n\t\t ax: {np.round(a[0],5)} [m/s2]'
    f'\n\t\t ay: {np.round(a[1],5)} [m/s2]'
    f'\n\t\t az: {np.round(a[2],5)} [m/s2]'
    f'\n\tThe maxmimu disturbance torques due to gust load in x,y,z'
    f'\n\t\t Mx, My, Mz: {np.round(M,2)}'
    f'\n\tThe maxmimum angular acceleration due to distrubance torque in x,y,z'
    f'\n\t\t alphax: {np.round(alpha[0],5)} [rad/s2]'
    f'\n\t\t alphay: {np.round(alpha[1],5)} [rad/s2]'
    f'\n\t\t alphaz: {np.round(alpha[2],5)} [rad/s2]'
)

