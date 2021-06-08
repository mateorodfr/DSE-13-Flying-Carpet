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

if section.shape == "square":
    dxy = 0.001
    Vymax, Vxmax = 1 , 1
    xyarr = np.arange(0,section.w * 0.5 +section.h * 0.5 + dxy, dxy)
    shearVy = (Vymax/section.Ix) * (section.h/2) * section.t_w * (xyarr < section.w/2) * xyarr + \
    (Vymax/section.Ix) * (section.h/2) * section.t_w * (xyarr >= section.w/2) * (section.w/2) + \
    (Vymax/section.Ix) * (section.t_h / 2) * 0.5*(xyarr - section.w/2) ** 2 * (xyarr > (section.w/2))

    #shearVx = (Vxmax/section.Iyy) * (section.w/2)*

    fig, ax = plt.subplots(1,1)
    #ax = ax.ravel()
    fig.suptitle('Shear stresses around the crossection')
    print(xyarr,shearVy)

    ax.plot(xyarr, shearVy)
    ax.title.set_text(r"Shear force in y along xy")
    ax.set_ylabel(r'tauy - Shear force in xy [N]')
    ax.set_xlabel(r'xy [m]')
    plt.show()