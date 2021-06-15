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
        if component == 'rotor':
            fig, ax = plt.subplots(3,4)
            ax = ax.ravel()
            fig.suptitle('Bridge internal loadings and deflections')

            ax[0].plot(z, np.flip(Vy_int))
            ax[0].title.set_text(r"Shear force in y along z")
            ax[0].set_ylabel(r'Vy - Shear force in y [N]')
            ax[0].set_xlabel(r'z [m]')

            ax[1].plot(z, np.flip(Mx_int))
            ax[1].title.set_text(r"Moment around x along z")
            ax[1].set_ylabel(r'Mx - Moment around x [Nm]')
            ax[1].set_xlabel(r'z [m]')

            ax[2].plot(z, np.flip(thetay))
            ax[2].title.set_text(r"Deflection angle in y along z")
            ax[2].set_ylabel(r'$\theta_y$ - Deflection angle in y [rad]')
            ax[2].set_xlabel(r'z [m]')

            ax[3].plot(z, np.flip(deflecty))
            ax[3].title.set_text(r"Deflection in y along z")
            ax[3].set_ylabel(r'vy - Delfection [m]')
            ax[3].set_xlabel(r'z [m]')

            ax[4].plot(z, np.flip(Vx_int))
            ax[4].title.set_text(r"Shear force in x along z")
            ax[4].set_ylabel(r'Vx - Shear force in x [N]')
            ax[4].set_xlabel(r'z [m]')

            ax[5].plot(z, np.flip(My_int))
            ax[5].title.set_text(r"Moment around y along z")
            ax[5].set_ylabel(r'My - Moment around y [Nm]')
            ax[5].set_xlabel(r'z [m]')

            ax[6].plot(z, np.flip(thetax))
            ax[6].title.set_text(r"Deflection angle in x along z")
            ax[6].set_ylabel(r'$\theta_x$ - Deflection angle in x [rad]')
            ax[6].set_xlabel(r'z [m]')

            ax[7].plot(z, np.flip(deflectx))
            ax[7].title.set_text(r"Deflection in x along z")
            ax[7].set_ylabel(r'vx - Delfection [m]')
            ax[7].set_xlabel(r'z [m]')

            ax[8].plot(z, np.flip(Nz_int))
            ax[8].title.set_text(r"Normal force in z along z")
            ax[8].set_ylabel(r'Nz - Normal force in z [N]')
            ax[8].set_xlabel(r'z [m]')

            ax[9].plot(z, np.flip(Tz_int))
            ax[9].title.set_text(r"Torque around z along z")
            ax[9].set_ylabel(r'Tz - Torque around z [Nm]')
            ax[9].set_xlabel(r'z [m]')

            ax[10].plot(z, np.flip(dtheta))
            ax[10].title.set_text(r"Rate of twist along z")
            ax[10].set_ylabel(r'$d\theta/dz$ - Rate of twist[rad/m]')
            ax[10].set_xlabel(r'z [m]')

            ax[11].plot(z, np.flip(thetaz))
            ax[11].title.set_text(r"Twist angle")
            ax[11].set_ylabel(r'$\theta_z$ - Twist [rad]')
            ax[11].set_xlabel(r'z [m]')

            plt.show()

        elif component == 'dick':
            fig, ax = plt.subplots(2,2)
            ax = ax.ravel()
            fig.suptitle('Rotor arm internal loadings')

            ax[0].plot(z, np.flip(thetay))
            ax[0].title.set_text(r"Deflection angle in y along z")
            ax[0].set_ylabel(r'$d\theta/dz$ - Rate of twist[rad/m]')
            ax[0].set_xlabel(r'z [m]')

            ax[1].plot(z, np.flip(deflecty))
            ax[1].title.set_text(r"Deflection in y along z")
            ax[1].set_ylabel(r'vy - Delfection [m]')
            ax[1].set_xlabel(r'z [m]')

            ax[2].plot(z, np.flip(thetax))
            ax[2].title.set_text(r"Deflection angle in x along z")
            ax[2].set_ylabel(r'$d\theta/dz$ - Rate of twist[rad/m]')
            ax[2].set_xlabel(r'z [m]')

            ax[3].plot(z, np.flip(deflectx))
            ax[3].title.set_text(r"Deflection in x along z")
            ax[3].set_ylabel(r'vx - Delfection [m]')
            ax[3].set_xlabel(r'z [m]')
            plt.show()


        elif component == 'cheese':
            fig, ax = plt.subplots(2,2)
            ax = ax.ravel()
            fig.suptitle('Rotor arm internal loadings')

            ax[0].plot(z, np.flip(Vy_int))
            ax[0].title.set_text(r"Shear force in y along z")
            ax[0].set_ylabel(r'Vy - Shear force in y [N]')
            ax[0].set_xlabel(r'z [m]')

            ax[1].plot(z, np.flip(Mx_int))
            ax[1].title.set_text(r"Moment around x along z")
            ax[1].set_ylabel(r'Mx - Moment around x [Nm]')
            ax[1].set_xlabel(r'z [m]')

            ax[2].plot(z, np.flip(Vx_int))
            ax[2].title.set_text(r"Shear force in x along z")
            ax[2].set_ylabel(r'Vx - Shear force in x [N]')
            ax[2].set_xlabel(r'z [m]')

            ax[3].plot(z, np.flip(My_int))
            ax[3].title.set_text(r"Moment around y along z")
            ax[3].set_ylabel(r'My - Moment around y [Nm]')
            ax[3].set_xlabel(r'z [m]')
            
            plt.show()
            print(np.min(Vy_int),np.max(Mx_int),np.min(Vx_int),np.max(My_int))



    return reactions, z, Vy_int, Mx_int, thetay, deflecty, Vx_int, My_int, thetax, deflectx, Nz_int, Tz_int, dtheta, thetaz

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

def generateCrossSections(part,shape,prop_lim,t_lim,L,rho,syield,tauyield,deflectmax,sf,n=100,dz=0.01):

    SF = sf
    sigma_good = []
    ms = []
    sigma_maxs = []
    progress = -1
    sections = []
    tau_maxs = []


    if shape == 'circle':
        r = np.arange(prop_lim[0]/n,prop_lim[0]+2*prop_lim[0]/n,prop_lim[0]/n)
        t = np.arange(0.001,t_lim[0]+t_lim[0]/n, t_lim[0]/n)
        for i in range(len(r)):
            for j in range(len(t)):
                section = pm.CrossSectionParameters(shape,[r[i]],[t[j]])
                m = section.A*L*rho
                w0,P,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz = getReactions(part,section,rho)
                reactions, z, Vy_int, Mx_int, thetay, deflecty, Vx_int, My_int, thetax, deflectx, Nz_int, Tz_int, dtheta, thetaz = internalLoading(dz,section,mc,L, Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0)
                sigmas, sigma_max, _ = normalStress(section, Nz_int, Mx_int, My_int,z)

                _, taumax, _, _ = internalShear(section,np.amax(np.abs(Vy_int)),np.amax(np.abs(Vx_int)),np.amax(np.abs(Tz_int)))

                if np.abs(sigma_max*SF) <= syield and taumax*SF <= tauyield and np.max(np.abs(deflectx)) <= deflectmax and np.max(np.abs(deflecty)) <= deflectmax:
                    sigma_good.append(sigmas)
                    ms.append(m)
                    sigma_maxs.append(sigma_max)
                    tau_maxs.append(taumax)
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
                    reactions, z, Vy_int, Mx_int, thetay, deflecty, Vx_int, My_int, thetax, deflectx, Nz_int, Tz_int, dtheta, thetaz = internalLoading(dz,section,mc,L, Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0)
                    sigmas, sigma_max, _ = normalStress(section, Nz_int, Mx_int, My_int,z)
                    _, taumax, _, _ = internalShear(section,np.amax(np.abs(Vy_int)),np.amax(np.abs(Vx_int)),np.amax(np.abs(Tz_int)))
                    if np.abs(sigma_max*SF) <= syield and taumax*SF <= tauyield and np.max(np.abs(deflectx)) <= deflectmax and np.max(np.abs(deflecty)) <= deflectmax:
                        sigma_good.append(sigmas)
                        ms.append(m)
                        sigma_maxs.append(sigma_max)
                        tau_maxs.append(taumax)
                        sections.append(section)
            if np.round((i/h[-1])*100) > progress:
                progress = np.round((i/h[-1])*100)
                print(f'Progress: {progress} %')






    return np.array(sigma_good),np.array(ms),np.array(sigma_maxs),sections, tau_maxs

def getReactions(part,section,rho):

    if part == 'bridge':
        w0 = rho * section.A * 9.80665
        P = 100*9.80665
        Fdx,Fdy,Fdz = 203,203,203
        Mdx,Mdy,Mdz = 56,56,56

    elif part == 'rotor':
        w0 = rho * section.A * 9.80665
        P = 0
        Fdx,Fdy,Fdz = 0,-1.1*2400*9.80665/4,0
        Mdx,Mdy,Mdz = 0,650,0
    return w0,P,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz

def internalShear(section,Vymax,Vxmax,Tmax,dxy=0.001,dtheta=0.0001,plot=False):


    if section.shape == "square":
        dxy = 0.001
        xyarr = np.arange(0,section.w * 0.5 +section.h * 0.5 + dxy, dxy)
        #tauy = ( ( (Vymax/section.Ix) * (section.h/2) * section.t_w * (xyarr < section.w/2) * xyarr + \
        #           (Vymax/section.Ix) * (section.h/2) * section.t_w * (xyarr >= section.w/2) * (section.w/2) )/section.t_w + \
        #((Vymax/section.Ix) * (section.t_h / 2) *( -0.5*(xyarr - section.w/2) ** 2  + section.h*(xyarr - section.w/2)/2 ) * (xyarr > (section.w/2)))/section.t_h)
        tauy = (Vymax/section.Ix)*section.t_w * (section.h/2) * xyarr *(xyarr <= section.w/2)/section.t_w + \
            ((Vymax/section.Ix)* section.t_w * (section.h/2) * (section.w/2) *(xyarr > section.w/2) + \
             ((Vymax/section.Ix) * section.t_h *   ((section.h/2) * (xyarr - section.w/2)*(xyarr > section.w/2)  - (1/2)*(xyarr-section.w/2)**2 *(xyarr > section.w/2)  )))/section.t_h 
             
             # (1/2) * (xyarr - section.w/2) ** 2 *(xyarr > section.w/2)))/section.t_h
       
        #taux = np.flip(((Vxmax/section.Iy) * (section.w/2) * section.t_h * (xyarr < section.h/2) *xyarr + \
        #((Vxmax/section.Iy) * (section.w/2) * section.t_h) * (xyarr >= section.h/2)  * section.h/2)/section.t_h + \
        #((Vxmax/section.Ix) * (section.t_w / 2) * (-0.5*(xyarr - section.h/2) ** 2 + section.w*(xyarr - section.h/2)/2 ) * (xyarr > (section.h/2)))/section.t_w)
        taux = np.flip((Vxmax/section.Iy)*section.t_h * (section.w/2) * xyarr *(xyarr <= section.h/2)/section.t_h + \
            ((Vxmax/section.Iy)* section.t_h * (section.w/2) * (section.h/2) *(xyarr > section.h/2) + \
             ((Vxmax/section.Iy) * section.t_w *   ((section.w/2) * (xyarr - section.h/2)*(xyarr > section.h/2)  - (1/2)*(xyarr-section.h/2)**2 *(xyarr > section.h/2)  )))/section.t_w )

        tauz = Tmax / (2*section.Am) * ((xyarr < (section.w /2)) * (1/section.t_w)   +  (xyarr >= (section.w /2)) * (1/section.t_h))

        #What I am doing here is sooooo sketch   
        tauy = np.concatenate((tauy,np.flip(tauy)[1:],-tauy[1:],np.flip(-tauy)[1:]))
        taux = np.concatenate((taux,np.flip(taux)[1:],-taux[1:],np.flip(-taux)[1:]))
        tauz = np.concatenate((tauz,tauz[1:],tauz[1:],tauz[1:]))
        xyarr = np.arange(0,section.w *2 +section.h * 2 +4* dxy, dxy)


        tau = taux + tauy + tauz

    elif section.shape == "circle":

        dtheta = 0.0001
        xyarr = np.arange(0, np.pi /2 + dtheta, dtheta)

        #tauy = ((Vymax/section.Ix) * section.t * section.r * ( - np.sin(xyarr / (2*np.pi*section.r)) * 1/(2*np.pi*section.r)))/section.t
        #tauy = (Vymax/section.Ix) * section.t * section.r * np.sin(xyarr)
        tauy = Vymax * np.sin(xyarr) / (np.pi * section.r * section.t)
        

        #taux = np.flip((Vxmax/section.Iy) * section.t * section.r * ( - np.sin(xyarr / (2*np.pi*section.r)) * 1/(2*np.pi*section.r)))/section.t
        #taux = np.flip((Vymax/section.Ix) * section.t * section.r * np.sin(xyarr))
        taux = np.flip(Vxmax * np.sin(xyarr) / (np.pi * section.r * section.t))

        tauz = (Tmax / (2*section.Am * section.t)) * np.ones(len(xyarr))

        tau = taux + tauy + tauz
        
    elif section.shape == "ibeam":
        raise NotImplementedError("Shear calculations for open sections not implemented")

    if plot:
        fig, ax = plt.subplots(2,2)
        ax = ax.ravel()
        fig.suptitle('Shear stresses around the crossection')

        ax[0].plot(xyarr, taux/1e6)
        ax[0].title.set_text(r"Shear stress caused by Vx")
        ax[0].set_ylabel(r'taux - Shear stress in xy [MPa]')
        ax[0].set_xlabel(r'xy [m]')

        ax[1].plot(xyarr, tauy/1e6)
        ax[1].title.set_text(r"Shear stress caused by Vy")
        ax[1].set_ylabel(r'tauy - Shear stress in xy [MPa]')
        ax[1].set_xlabel(r'xy [m]')

        ax[2].plot(xyarr, tauz/1e6)
        ax[2].title.set_text(r"Shear stress caused by T")
        ax[2].set_ylabel(r'tauz - Shear stress in xy [MPa]')
        ax[2].set_xlabel(r'xy [m]')

        ax[3].plot(xyarr, tau/1e6)
        ax[3].title.set_text(r"Total shear stress")
        ax[3].set_ylabel(r'tau - Shear stress in xy [MPa]')
        ax[3].set_xlabel(r'xy [m]')

        plt.show()

    return tau, np.amax(tau), [taux,tauy,tauz], xyarr




# ALU 7075-T651   E = 71.7 GPA G = 27.5 GPa ystress = 500MPA  shearstress  = 331 MPa rho = 2700 Cost = 1.8 $/kg

# ALU 2024-T4     E = 73.1 GPa G = 28 GPa  ystress = 324 MPa shearstress = 283 MPa rho = 2770 Cost = 2 $/kg

# TI 6Al4V        E = 114 GPA G = 40 GPA  ystress = 1100 MPA  shearstress = 550 #calc    rho = 4.5 cost = 16.25

# Steel, high strength 4340  E = 210 GPA G = 76 ystress = 1240  MPA shearstress = 620 #calc rho = 7.8 cost = 0.25

# Steel, mild 1020 E = 210 GPA G = 76 GPS ystress = 200 shearstress = 100 rho = 7.8 cost = 0.5

rho_bridge = 2700
w_bridge = 0.1
h_bridge = 0.1
thicc_w = 0.001
thicc_h = 0.001
L_bridge = 1.72

# section = pm.CrossSectionParameters('square', [h_bridge,w_bridge],[thicc_h,thicc_w])

E = 71.7e9
G = 27.5e9
Money = 1.8
dz = 0.01


plotInternal = True


sigma_yield = 500e6
tau_yield = 331e6
deflectmax = 0.02
component = 'rotor'
r_lim = [0.2]
t_lim = [0.01]
hw_lim = [0.3,0.1]
SF = 1.5
n = 10



# w0,P,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz = getReactions(component, section,rho_bridge)
# reactions,z, Vy_int, Mx_int, thetay , deflecty, Vx_int, My_int, thetax, deflectx, Nz_int, Tz_int, dtheta, thetaz = internalLoading(dz,section,mc,L_bridge, Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0,plotInternal)
# sigmas, sigma_max, z_max = normalStress(section, Nz_int, Mx_int, My_int,z)
# tau, taumax, tauxyz,xyarr = internalShear(section, np.amax(np.abs(Vy_int)),np.amax(np.abs(Vx_int)),np.amax(np.abs(Tz_int)),plot=True)
sigma_good,ms,sigma_maxs,sections, tau_maxs = generateCrossSections(component,'square',hw_lim,t_lim,L_bridge,rho_bridge,sigma_yield,tau_yield,deflectmax,SF,n)



idx = np.where(ms == np.amin(ms))[0][0]

section_best = sections[idx]
section_best.plotNormalStress(sigma_good[idx],tau_maxs[idx],ms[idx],sigma_yield)

w0,P,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz = getReactions('rotor',section_best,rho_bridge)
reactions, z, Vy_int, Mx_int, thetay, deflecty, Vx_int, My_int, thetax, deflectx, Nz_int, Tz_int, dtheta, thetaz=internalLoading(dz,section_best,mc,L_bridge,Fdx,Fdy,Fdz,Mdx,Mdy,Mdz,P,w0,plot=True)

tau, maxtau, [taux,tauy,tauz], xyarr  = internalShear(section_best,-1.1*2400*9.80665/4,0,0,dxy=0.001,dtheta=0.0001,plot=True)



print(section_best.t_h, section_best.Ix)
print(len(sections))