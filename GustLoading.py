import numpy as np
import Parameters as pm
import matplotlib.pyplot as plt

def dragForce(V,S,Cd,rho):
    return [0.5*rho*S[0]*Cd*V**2,0.5*rho*S[1]*Cd*V**2,0.5*rho*S[2]*Cd*V**2]

def gustDisturbance(S,Uds,t,drag,dt,Cd,rho):

    ts = np.arange(0,2*t+dt,dt)
    U = Uds/2*(1-np.cos((np.pi*ts)/(ts[-1]/2)))
    D = drag(U,S,Cd,rho)

    return np.transpose(np.array(D)),np.array(U),np.array(ts)

def gustEnvelope(Hmax,hmax,Umax,concept,plot=False):

    H = np.arange(10,Hmax,Hmax/100)
    Fgz = 1-hmax/76200
    R1 = concept.Mtot_concept/(concept.Mtot_concept-500)
    R2 = 1
    Fgm = np.sqrt(R2*np.tan(np.pi*R1/4))
    Fg = 0.5*(Fgz+Fgm)

    Us = []
    Ss = []
    Udss = []
    for i in range(1,len(H)):
        Uds = Umax*Fg*(H[i]/107)**(1/6)
        s = np.arange(0,2*H[i],1)
        U = Uds/2*(1-np.cos(np.pi*s/(s[-1]/2)))
        Us.append(U*0.3048)
        Udss.append(Uds*0.3048)
        Ss.append(s*0.3048)
        if plot:
            plt.plot(Ss[i-1],Us[i-1])
    if plot:
        plt.title('Gust Envelope')
        plt.xlabel('Time t [s]')
        plt.ylabel(r'Gust intensity U [ms$^{-1}$]')
        plt.grid()
        plt.savefig('figures/gustenvelope.png')
        plt.show()

    maxUds = np.max(Udss)
    return maxUds

def gustMax(Ui,zi,hmax,z0,dz):
    dz = 0.5 #height step [m]
    z = np.arange(zi,zf+dz,dz) #operational height range [m]
    Uf = Ui*(np.log(z/z0)/np.log(zi/z0)) #Gust speed at given height [m/s]
    return Uf[-1],z,Uf #Max gust load

def getGustData(concept,tsim,Ui,zi,hmax,z0,Hmax,Cd=1.2,rho=1.225,dz=0.5,dtsim=0.1,plotGust=False,plotDisturbance=False,isPrint=False,saveText=False):

    #Define Masses
    m = concept.Mtot_concept-concept.Mpay_concept #curent final iterated mass

    #Cabin dimensions
    L_cabin = concept.cabin.L_cabin #Cabin length
    W_cabin = concept.cabin.W_cabin #Cabin width
    H_cabin = concept.cabin.H_cabin #Cabin height
    Dim_cabin = np.array([L_cabin,W_cabin,H_cabin]) #Array to store all cabin dimensions
    d = (1/6)*Dim_cabin #Moment arm in each direction

    Syz = concept.cabin.S_cabin[2] #yz plane
    Sxz = concept.cabin.S_cabin[1] #xz plane
    Sxy = concept.cabin.S_cabin[0] #xy plane
    S = np.array([Syz,Sxz,Sxy])

    #Mass Moment of Inertia
    I = concept.MMOI() #Mass moment of inertia stored n array xx,yy,zz


    """Max Gust Speed from Certification Model"""

    #Gust Certification parameters
    Hmt = Hmax
    Umax, z, Uf = gustMax(Ui,zi,hmax,z0,dz)
    Udsmax = gustEnvelope(Hmt,hmax,Umax,concept,plotGust)
    print(Udsmax)
    """
    Compute characteristics

    Reasoning for Trans. Accelerations
        The max gust speed is the same in all directions
        The gust speed taken is the max vertical gust speed.
        The max lateral gust speed is equal to vertical gust speed
        Gust drag force is distributed around cg such that no torque is produced

    Reasoning for Rot. Accelerations
        Assuming the same magnitude of gust but distributed in a triangle.
        This will cause a translational acceleration equal to ax,ay,az calculated above.
        The distributed load will cause a torque.
        The load is simplified as a point load applied at 1/3 of the load.
        Assuming the load goes from 0 to w0 over the cabin
        This means the load is applied at 1/3 Lcabin, Wcabin, Hcabin
    """

    #Gust force over time
    dt = dtsim
    t = tsim
    Dgust,Ugust,tgust = gustDisturbance(S,Udsmax,t,dragForce,dt,Cd,rho)

    #Translational accelerations gust point load applied at cg
    a = Dgust/m

    #Rotational accelerations due to distributed gust load.
    M  = Dgust*d #Torque
    alpha = M/I #Angular accelerations


    """ Plot gusts"""

    if plotGust:
        plt.plot(z,Uf)
        plt.title('Gust speed at varying altitudes')
        plt.xlabel('Altitude: h [m]')
        plt.ylabel(r'Wind speed: U [ms$^{-1}$]')
        plt.grid()
        plt.savefig('figures/gustload.png')
        plt.show()
        # plt.savefig("figures/gustload")


    """Plot Disturbance Force"""

    if plotDisturbance:

        plt.plot(tgust,Dgust)
        plt.xlabel('t time [s]')
        plt.ylabel('Gust force in [N]')
        plt.title('Disturbance')
        plt.show()


    """Print Characteristics"""

    if isPrint:
        print(
            f'The gust load analysis yields: '
            f'\n\tThe moment of inertias in x,y,z'
            f'\n\t\t Ixx: {np.round(I[0])} [kgm2]'
            f'\n\t\t Iyy: {np.round(I[1])} [kgm2]'
            f'\n\t\t Izz: {np.round(I[2])} [kgm2]'
            f'\n\tThe maximum gust speed in x,y,z:'
            f'\n\t\t Umax: {np.round(Umax,2)} [m/s]'
            f'\n\tThe maximum distrubance force due to gust load in x,y,z:'
            f'\n\t\t Fx: {np.round(np.max(Dgust[:,0]),2)} [N]'
            f'\n\t\t Fy: {np.round(np.max(Dgust[:,1]),2)} [N]'
            f'\n\t\t Fz: {np.round(np.max(Dgust[:,2]),2)} [N]'
            f'\n\tThe maximum translational accelerations due to disturbance in x,y,z'
            f'\n\t\t ax: {np.round(np.max(a[:,0]),5)} [m/s2]'
            f'\n\t\t ay: {np.round(np.max(a[:,1]),5)} [m/s2]'
            f'\n\t\t az: {np.round(np.max(a[:,2]),5)} [m/s2]'
            f'\n\tThe maxmimu disturbance torques due to gust load in x,y,z'
            f'\n\t\t Mx: {np.round(np.max(M[:,0]),2)}'
            f'\n\t\t My: {np.round(np.max(M[:,1]),2)}'
            f'\n\t\t Mz: {np.round(np.max(M[:,2]),2)}'
            f'\n\tThe maxmimum angular acceleration due to distrubance torque in x,y,z'
            f'\n\t\t alphax: {np.round(np.max(alpha[:,0]),5)} [rad/s2]'
            f'\n\t\t alphay: {np.round(np.max(alpha[:,1]),5)} [rad/s2]'
            f'\n\t\t alphaz: {np.round(np.max(alpha[:,2]),5)} [rad/s2]'
        )

    if saveText:

        saveArrays(r,names,store)

    return Dgust,M,Ugust,tgust,a,alpha

def saveArrays(root,names,arr):

    for i in range(len(arr)):
        np.savetxt(root+names[i],arr[i])


""" Import General Parameters"""

#Import subsystem parameters
concept = pm.ConceptParameters(0)


#all these variable sneed to be manually changed to perform a case specific analysis
Cd = 1.05 #CD of squre at reynolds: 6e5
rho = concept.physics.rho0 #[kg/m3] density at sea level
tsim = 50
dtsim = 0.01


"""Define Module Specific parameters"""

#Gust load constants
Ui = 5 #wind speed at sensor height in [m/s]
z0 = 3 #surface roughness constant -> 3m for urban centers [m]
zi = 10 #Anemometer height when measuring wind spped of ground [m]
zf = 400 #Maximum operational height [m]
dz = 0.01
Hmax = 1000


""" Define Output Booleans"""

plotGust = False #if this is true the program will plot the gust speed as a fucntion of height and the gust envelope
plotDisturbance = False
isPrint = False #if this is true the program will print all the characteristics
saveText = True #Save arrays to text


r = r'Data/Gust/'
names = [r'Dgust.txt',r'Mgust.txt',r'Ugust.txt',r'tgust.txt',r'agust.txt',r'alphagust.txt']
Dgust,M,Ugust,tgust,a,alpha = getGustData(concept,tsim,Ui,zi,zf,z0,Hmax,Cd,rho,dz,dtsim,plotGust,plotDisturbance,isPrint)
store = [Dgust,M,Ugust,tgust,a,alpha]
saveArrays(r,names,store)
print(np.max(Dgust[:,0]))

