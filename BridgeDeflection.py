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

def normalStress(section, Nz, Mx, My,z):
    sigma_max = 0
    sigma_dist_max = []
    for i in range(len(Nz)):
        sigma_i = np.ones(len(section.contour))*(Nz[i]/section.A) \
                  + (My[i]/section.Iy) * section.contour[:,0]     \
                  + (Mx[i]/section.Ix) * section.contour[:,1]

        if np.max(np.abs(sigma_i)) > np.abs(sigma_max):
            sigma_low = np.amin(sigma_i)
            sigma_high = np.amax(sigma_i)
            if np.abs(sigma_low) > sigma_high:
                sigma_max = sigma_low
            else:
                sigma_max = sigma_high
            sigma_dist_max = sigma_i
            z_max = z[i]
    return sigma_dist_max, sigma_max, z_max

def generateCrossSections(part,shape,prop_lim,t_lim,L,rho,syield,sf,n=100,dz=0.01):

    SF = sf
    sigma_good = []
    ms = []
    sigma_maxs = []
    progress = -1
    sections = []


    if shape == 'circle':
        r = np.arange(prop_lim[0]/n,prop_lim[0]+2*prop_lim[0]/n,prop_lim[0]/n)
        t = np.arange(0.001,t_lim[0]+t_lim[0]/n, t_lim[0]/n)
        for i in range(len(r)):
            for j in range(len(t)):
                section = pm.CrossSectionParameters(shape,[r[i]],[t[j]])
                m = section.A*L*rho
                w0,P,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz = getReactions(part,section,rho)
                _,z, _, Mx_int, _, _, _, My_int, _, _, Nz_int, _, _, _ = internalLoading(dz,section,mc,L, Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0)
                sigmas, sigma_max, _ = normalStress(section, Nz_int, Mx_int, My_int,z)

                if np.abs(sigma_max*SF) <= syield:
                    sigma_good.append(sigmas)
                    ms.append(m)
                    sigma_maxs.append(sigma_max)
                    sections.append(section)


            if np.round((i/len(r))*100) > progress:
                progress = np.round((i/len(r))*100)
                print(f'Progress: {progress} %')

    elif shape == 'square':
        h = np.arange(prop_lim[0]/n,prop_lim[0]+2*prop_lim[0]/n,prop_lim[0]/n)
        w = np.arange(prop_lim[1]/n,prop_lim[1]+2*prop_lim[1]/n,prop_lim[1]/n)
        t = np.arange(0.001,t_lim[0]+2*t_lim[0]/(n/3),t_lim[0]/(n/3))

        for i in h:
            for j in w:
                for k in t:
                    section = pm.CrossSectionParameters(shape,[i,j],[k,k])
                    m = section.A*rho*L
                    w0,P,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz = getReactions(part,section,rho)
                    _,z, _, Mx_int, _, _, _, My_int, _, _, Nz_int, _, _, _ = internalLoading(dz,section,mc,L, Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0)
                    sigmas, sigma_max, _ = normalStress(section, Nz_int, Mx_int, My_int,z)
                    if np.abs(sigma_max*SF) <= syield:
                        sigma_good.append(sigmas)
                        ms.append(m)
                        sigma_maxs.append(sigma_max)
                        sections.append(section)
            if np.round((i/h[-1])*100) > progress:
                progress = np.round((i/h[-1])*100)
                print(f'Progress: {progress} %')






    return np.array(sigma_good),np.array(ms),np.array(sigma_maxs),sections

def getReactions(part,section,rho):

    if part == 'bridge':
        w0 = rho * section.A * 9.80665
        P = 100*9.80665
        Fdx,Fdy,Fdz = 203,203,203
        Mdx,Mdy,Mdz = 56,56,56

    elif part == 'rotor':
        w0 = rho * section.A * 9.80665
        P = 0
        Fdx,Fdy,Fdz = 0,2400*9.80665/4,0
        Mdx,Mdy,Mdz = 0,750,0
    return w0,P,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz

rho_bridge = 2700
w_bridge = 0.5
h_bridge = 0.1
thicc_w = 0.001
thicc_h = 0.001
L_bridge = 2.5

section = pm.CrossSectionParameters('square', [h_bridge,w_bridge],[thicc_h,thicc_w])

E = 69e9
G = 25.5e9
dz = 0.01

plotInternal = False


sigma_yield = 276e6
component = 'rotor'
r_lim = [1]
t_lim = [0.005]
hw_lim = [0.05,0.1]
SF = 1.5
n = 100


w0,P,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz = getReactions(component, section,rho_bridge)
reactions,z, Vy_int, Mx_int, thetay, deflecty, Vx_int, My_int, thetax, deflectx, Nz_int, Tz_int, dtheta, thetaz = internalLoading(dz,section,mc,L_bridge, Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0,plotInternal)
sigmas, sigma_max, z_max = normalStress(section, Nz_int, Mx_int, My_int,z)

sigma_good,ms,sigma_maxs,sections = generateCrossSections(component,'circle',r_lim,t_lim,L_bridge,rho_bridge,sigma_yield,SF,n)
idx = np.where(ms == np.amin(ms))[0][0]

section_best = sections[idx]
section_best.plotNormalStress(sigma_good[idx],ms[idx],sigma_yield)



# np.savetxt('Data/Stress/sigmas_good.txt', sigma_good)
# np.savetxt('Data/Stress/ms.txt', ms)
# np.savetxt('Data/Stress/sigmas_maxs.txt', sigma_maxs)
# np.savetxt('Data/Stress/props.txt', properties)












