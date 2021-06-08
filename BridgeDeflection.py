from Macaulay import macaulay as mc
import numpy as np
import matplotlib.pyplot as plt
import Parameters as pm

def internalLoading(dz,section,mc,L_beam,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0,plot=False):
    """Returns internal loading diagrams in x,y,z for given loading"""

    if section.shape == 'square':
        st_integral = 2*section.h/section.t_h + 2*section.w/section.t_w
    elif section.shape == 'circle':
        st_integral = (2*np.pi*section.r)/section.t
    elif section.shape == 'ibeam':
        st_integral = 2*section.w*section.t_w + section.h*section.t_h


    z = np.arange(0,L_beam+dz,dz)

    Ry = Fdy - 3 * P - w0 * L_beam
    Rx = Fdx
    Rz = Fdz

    Mrx = Ry*L_beam+P*(L_beam+L_beam/2)+0.5*w0*L_beam*L_beam + Mdx
    Mry = Rx*L_beam - Mdy
    Mrz = Mdz

    reactions = [Ry,Rx,Rz,Mrx,Mry,Mrz]

    Vy_int = -Ry*mc(z,0,0) -P*mc(z,0,0) - w0*mc(z,0,1) - P*mc(z,L_beam/2,0) - P*mc(z,L_beam,0) + Fdy*mc(z,L_beam,0)
    Mx_int = -Ry*mc(z,0,1) - P*mc(z,0,1) - 0.5*w0*mc(z,0,2) - P*mc(z,L_beam/2,1) + Mrx*mc(z,0,0)
    thetay = (-1/(E*section.Ix)) * ( -0.5*Ry*mc(z,0,2) -0.5*P*mc(z,0,2) -(1/6)*w0*mc(z,0,3) -0.5*P*mc(z,L_beam/2,2) + Mrx*mc(z,0,1)    )
    deflecty = (-1/(E*section.Ix)) * ( (-1/6) * Ry * mc(z,0,3) + (-1/6) * P * mc(z,0,3) + (-1/24) * w0 * mc(z,0,4) + (1/2) * Mrx * mc(z,0,2) + (-1/6) * P *mc(z,L_beam,3) )

    Vx_int = -Rx*mc(z,0,0) + Fdx*mc(z,L_beam,0)
    My_int =  -Rx * mc(z,0,1) + Mry * mc(z,0,0)
    thetax = (-1/(E*section.Iy))* ((-1/2) * Rx * mc(z,0,2) + Mry * mc(z,0,1))
    deflectx = (-1/(E*section.Iy)) * ((-1/6)*Rx*mc(z,0,3) + (1/2)*Mry*mc(z,0,2))

    Nz_int = - Rz * mc(z,0,0) + Fdz*mc(z,L_beam,0)
    Tz_int = Mrz * mc(z,0,0) - Mdz*mc(z,L_beam,0)

    if section.shape != 'ibeam':
        dtheta = st_integral*Tz_int/(4 * section.Am**2 * G)
        thetaz = np.cumsum(dtheta*dz)
    else:
        raise NotImplementedError('Twist for open sections needs to be implemented')
        dtheta, thetaz = None, None
    if plot:
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

        ax[10].plot(z, dtheta)
        ax[10].title.set_text(r"Rate of twist along z")
        ax[10].set_ylabel(r'$d\theta/dz$ - Rate of twist[rad/m]')
        ax[10].set_xlabel(r'z [m]')

        ax[11].plot(z, thetaz)
        ax[11].title.set_text(r"Twist angle")
        ax[11].set_ylabel(r'$\theta_z$ - Twist [rad]')
        ax[11].set_xlabel(r'z [m]')

        plt.show()

    return reactions,z, Vy_int, Mx_int, thetay, deflecty, Vx_int, My_int, thetax, deflectx, Nz_int, Tz_int, dtheta, thetaz


rho_bridge = 2700
w_bridge = 0.5
h_bridge = 0.1
thicc_w = 0.001
thicc_h = 0.001
section = pm.CrossSectionParameters('square',[h_bridge,w_bridge],[thicc_h,thicc_w])
sectionn = pm.CrossSectionParameters('circle',[h_bridge],[thicc_h])

Ix = section.Ix
Iy = section.Iy
Jz = section.Jz


E = 69e9
G = 25.5e9
L_bridge = 2.5

isBridge = False
isRotor = True

if isBridge:
    w0 = rho_bridge * section.A * 9.80665
    P = 100*9.80665
    Fdx,Fdy,Fdz = 203,203,203
    Mdx,Mdy,Mdz = 56,56,56
elif isRotor:
    w0 = rho_bridge * section.A * 9.80665
    P = 0
    Fdx,Fdy,Fdz = 0,2400*9.80665/4,0
    Mdx,Mdy,Mdz = 0,750,0

#Set up reactions
dz = 0.01
plotInternal = False
reactions,z, Vy_int, Mx_int, thetay, deflecty, Vx_int, My_int, thetax, deflectx, Nz_int, Tz_int, dtheta, thetaz = internalLoading(dz,section,mc,L_bridge, Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0,plotInternal)


def normalStress(Nz,My,Mx,section):

    sigma = Nz/section.A + (Mx/section.Ix)*section.contour[:,1] + (My/section.Iy)*section.contour[:,0]
    print(len(sigma))
normalStress(Nz_int,My_int,Mx_int,section)