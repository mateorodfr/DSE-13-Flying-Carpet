import Parameters as pm
import numpy as np
from BridgeDeflection import section #,Nz_int,Mx_int,My_int
import matplotlib.pyplot as plt

Ixx = section.Ix
Iyy = section.Iy

ys = [-section.h/2,section.h/2,section.h/2,-section.h/2]
xs = [-section.w/2,section.w/2,-section.w/2,section.w/2]





# sigmas = []
# maxsigma = []
# for i in range(len(ys)):
#     sigma_z = Nz_int/section.A + Mx_int/Ixx * ys[i] + (My_int/Iyy) * xs[i]
#     sigmas.append(sigma_z)

#     if np.sum(np.abs(sigma_z)) > np.sum(np.abs(maxsigma)):
#         maxsigma =  np.sum(np.abs(sigma_z))
#         sigma_final = sigma_z



taus = []
maxtau = []

dxy = 0.001
Vymax, Vxmax , Tmax = 1 , 1, 1

if section.shape == "square":
    xyarr = np.arange(0,section.w * 0.5 +section.h * 0.5 + dxy, dxy)
    tauy = ((Vymax/section.Ix) * (section.h/2) * section.t_w * (xyarr < section.w/2) * xyarr + \
    (Vymax/section.Ix) * (section.h/2) * section.t_w * (xyarr >= section.w/2) * (section.w/2))/section.t_w + \
    ((Vymax/section.Ix) * (section.t_h / 2) *( -0.5*(xyarr - section.w/2) ** 2  +section.h*(xyarr - section.w/2)/2 ) * (xyarr > (section.w/2)))/section.t_h

    taux = np.flip(((Vxmax/section.Iy) * (section.w/2) * section.t_h * (xyarr < section.h/2) *xyarr + \
    ((Vxmax/section.Iy) * (section.w/2) * section.t_h) * (xyarr >= section.h/2)  * section.h/2)/section.t_h + \
    ((Vxmax/section.Ix) * (section.t_w / 2) * (-0.5*(xyarr - section.h/2) ** 2 + section.w*(xyarr - section.h/2)/2 ) * (xyarr > (section.h/2)))/section.t_w)

    tauz = Tmax / (2*section.Am) * ( (xyarr < (section.w /2)) * (1/section.t_w)   +  (xyarr >= (section.w /2)) * (1/section.t_h)    )

    tau = taux + tauy + tauz
    #shearVx = (Vxmax/section.Iyy) * (section.w/2)*

if section.shape == "circle"
    xyarr = np.arange(0,section.r * )


fig, ax = plt.subplots(2,2)
ax = ax.ravel()
fig.suptitle('Shear stresses around the crossection')

ax[0].plot(xyarr, taux)
ax[0].title.set_text(r"Shear stress caused by Vx")
ax[0].set_ylabel(r'tauy - Shear flow in xy [N]')
ax[0].set_xlabel(r'xy [m]')

ax[1].plot(xyarr, tauy)
ax[1].title.set_text(r"Shear stress caused by Vy")
ax[1].set_ylabel(r'tauy - Shear flow in xy [N]')
ax[1].set_xlabel(r'xy [m]')

ax[2].plot(xyarr, tauz)
ax[2].title.set_text(r"Shear stress caused by T")
ax[2].set_ylabel(r'tauy - Shear flow in xy [N]')
ax[2].set_xlabel(r'xy [m]')

ax[3].plot(xyarr, tau)
ax[3].title.set_text(r"Total shear stress")
ax[3].set_ylabel(r'tauy - Shear flow in xy [N]')1
ax[3].set_xlabel(r'xy [m]')

plt.show()