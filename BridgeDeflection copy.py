from Macaulay import macaulay as mc
import numpy as np
import matplotlib.pyplot as plt




rho_bridge = 2700
w_bridge = 0.5
h_bridge = 0.1
thicc_w = 0.001
thicc_h = 0.001
thicc_avg = (thicc_h+thicc_w)/2


A_bridge = 2*w_bridge * thicc_w + 2* h_bridge * thicc_h
A_m = w_bridge*h_bridge
st_integral = 2*h_bridge/thicc_h + 2*w_bridge/thicc_w

Ixx = (1/3) * h_bridge**2 * w_bridge * thicc_avg
Iyy = (1/3) * w_bridge**2 * h_bridge * thicc_avg
Jz = (w_bridge*h_bridge*thicc_avg/3)*(h_bridge+w_bridge)



E = 69e9
G = 25.5e9
L_bridge = 2.5


w0 = rho_bridge * A_bridge * 9.80665

P = 100*9.80665
Fdx,Fdy,Fdz = 203,203,203
Mdx,Mdy,Mdz = 56,56,56

dz = 0.01
z = np.arange(0,L_bridge+dz,dz)

def internalloadings(E,Ixx,Iyy,z,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0,L_bridge):
    #Set up reactions
    Ry = Fdy -3*P-w0*L_bridge
    Rx = Fdx
    Rz = Fdz

    Mrx = Ry*L_bridge+P*(L_bridge+L_bridge/2)+0.5*w0*L_bridge*L_bridge + Mdx
    Mry = Rx*L_bridge - Mdy
    Mrz = Mdz

    #Axis Y
    Vy_int = -Ry*mc(z,0,0) -P*mc(z,0,0) - w0*mc(z,0,1) - P*mc(z,L_bridge/2,0) - P*mc(z,L_bridge,0) + Fdy*mc(z,L_bridge,0)
    Mx_int = -Ry*mc(z,0,1) - P*mc(z,0,1) - 0.5*w0*mc(z,0,2) - P*mc(z,L_bridge/2,1) + Mrx*mc(z,0,0)
    thetay = (-1/(E*Ixx)) * ( -0.5*Ry*mc(z,0,2) -0.5*P*mc(z,0,2) -(1/6)*w0*mc(z,0,3) -0.5*P*mc(z,L_bridge/2,2) + Mrx*mc(z,0,1)    )
    deflecty = (-1/(E*Ixx)) * ( (-1/6) * Ry * mc(z,0,3) + (-1/6) * P * mc(z,0,3) + (-1/24) * w0 * mc(z,0,4) + (1/2) * Mrx * mc(z,0,2) + (-1/6) * P *mc(z,L_bridge,3) )

    #Axis X
    Vx_int = -Rx*mc(z,0,0) + Fdx*mc(z,L_bridge,0)
    My_int =  -Rx * mc(z,0,1) + Mry * mc(z,0,0)
    thetax = (-1/(E*Iyy))* ((-1/2) * Rx * mc(z,0,2) + Mry * mc(z,0,1))
    deflectx = (-1/(E*Iyy)) * ((-1/6)*Rx*mc(z,0,3) + (1/2)*Mry*mc(z,0,2))

    #Axis Z
    Nz_int = - Rz * mc(z,0,0) + Fdz*mc(z,L_bridge,0)
    Tz_int = Mrz * mc(z,0,0) - Mdz*mc(z,L_bridge,0)
    dtheta = st_integral*Tz_int/(4 * A_m**2 * G)
    thetaz = np.cumsum(dtheta*dz)
    return Vy_int , Mx_int , thetay , deflecty , Vx_int , My_int , thetax , deflectx , Nz_int , Tz_int , dtheta , thetaz



Vy_int , Mx_int , thetay , deflecty , Vx_int , My_int , thetax , deflectx , Nz_int , Tz_int , dtheta , thetaz = internalloadings(E,Ixx,Iyy,z,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0,L_bridge)



ys = [-h_bridge/2,h_bridge/2,h_bridge/2,-h_bridge/2]
xs = [-w_bridge/2,w_bridge/2,-w_bridge/2,w_bridge]
sigmas = []
maxsigma = []
for i in range(len(ys)):
    sigma_z = Nz_int/A_bridge + Mx_int/Ixx * ys[i] + (My_int/Iyy) * xs[i]
    sigmas.append(sigma_z)
    if np.sum(np.abs(sigma_z)) > np.sum(np.abs(maxsigma)):
        maxsigma =  np.sum(np.abs(sigma_z))
        sigma_final = sigma_z

plt.plot(z,sigma_final)
plt.show()

plotInternal = True
if plotInternal:
    fig, ax = plt.subplots(3,4)
    ax = ax.ravel()
    fig.suptitle('Drawbridge deflection in Y')

    ax[0].plot(z, Vy_int)
    ax[0].title.set_text(r"Shear force in y along z")
    ax[0].set_ylabel(r'Vy - Shear force in y [N]')
    ax[0].set_xlabel(r'z [m]')

    ax[1].plot(z, Mx_int)
    ax[1].title.set_text(r"Moment around x along z")
    ax[1].set_ylabel(r'Mx - Moment around x [Nm]')
    ax[1].set_xlabel(r'z [m]')

    ax[2].plot(z, thetay)
    ax[2].title.set_text(r"Deflection angle in y along z")
    ax[2].set_ylabel(r'$\theta_y$ - Deflection angle in y [rad]')
    ax[2].set_xlabel(r'z [m]')

    ax[3].plot(z, deflecty)
    ax[3].title.set_text(r"Deflection in y along z")
    ax[3].set_ylabel(r'vy - Delfection [m]')
    ax[3].set_xlabel(r'z [m]')

    ax[4].plot(z, Vx_int)
    ax[4].title.set_text(r"Shear force in x along z")
    ax[4].set_ylabel(r'Vx - Shear force in x [N]')
    ax[4].set_xlabel(r'z [m]')

    ax[5].plot(z, My_int)
    ax[5].title.set_text(r"Moment around y along z")
    ax[5].set_ylabel(r'My - Moment around y [Nm]')
    ax[5].set_xlabel(r'z [m]')

    ax[6].plot(z, thetax)
    ax[6].title.set_text(r"Deflection angle in x along z")
    ax[6].set_ylabel(r'$\theta_x$ - Deflection angle in x [rad]')
    ax[6].set_xlabel(r'z [m]')

    ax[7].plot(z, deflectx)
    ax[7].title.set_text(r"Deflection in x along z")
    ax[7].set_ylabel(r'vx - Delfection [m]')
    ax[7].set_xlabel(r'z [m]')

    ax[8].plot(z, Nz_int)
    ax[8].title.set_text(r"Normal force in z along z")
    ax[8].set_ylabel(r'Nz - Normal force in z [N]')
    ax[8].set_xlabel(r'z [m]')

    ax[9].plot(z, Tz_int)
    ax[9].title.set_text(r"Torque around z along z")
    ax[9].set_ylabel(r'Tz - Torque around z [Nm]')
    ax[9].set_xlabel(r'z [m]')

    ax[10].plot(z, thetax)
    ax[10].title.set_text(r"Rate of twist along z")
    ax[10].set_ylabel(r'$d\theta/dz$ - Rate of twist[rad/m]')
    ax[10].set_xlabel(r'z [m]')

    ax[11].plot(z, deflectx)
    ax[11].title.set_text(r"Twist angle")
    ax[11].set_ylabel(r'$\theta_z$ - Twist [rad]')
    ax[11].set_xlabel(r'z [m]')

    plt.show()