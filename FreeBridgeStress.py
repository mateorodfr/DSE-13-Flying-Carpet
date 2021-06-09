import numpy as np
import matplotlib.pyplot as plt
import Parameters as pm
from mpl_toolkits.mplot3d import Axes3D

rho_bridge = 2700
E = 69e9
G = 25.5e9
L = 2.5
g=9.81

w_bridge = 0.5
h_bridge = 0.5
thicc_top = 0.001
thicc_bot = 0.001
thicc_side = 0.001

Volume= (2*thicc_top* w_bridge +2* thicc_side * h_bridge) *L
Total_weight= Volume * rho_bridge


props= [h_bridge,w_bridge]
ts=[thicc_side,thicc_top]

shape= pm.CrossSectionParameters('square',props,ts)

Ixx = shape.Ix


Fpay= 100*g
W_bridge= Total_weight/L *g


z= 0 #np.linspace(0,L,50)
count=0
dz= L/50

dlist=[]
Mlist=[]
vlist=[]
A= -Fpay/2*61/36*L**2 - W_bridge/6*L**3
B=-Fpay*307/6/216 *L**3 -W_bridge/24*L**4 - A*L
#v= -1/E/Ixx *(Fpay/6*z**3 + W_bridge/24 *z**4  - Fpay/2*L**2*z -W_bridge/6*L**3*z + Fpay/3*L**3 + W_bridge/8*L**4)

while z <= L:

    if z < L/3:
        deflection= -1/E/Ixx *(Fpay/6*z**3 + W_bridge/24 *z**4 + A*z +B)
        M = Fpay * z + W_bridge / 2 * z ** 2
        z += dz
        v= Fpay + W_bridge* z
        dlist.append(deflection)
        Mlist.append(M)
        vlist.append(v)


    elif z < L/2:
        deflection = -1 / E / Ixx * (Fpay / 6 * (z ** 3 + (z- 1/3*L)**3) + W_bridge / 24 * z ** 4 +A*z + B)
        M= Fpay*z +Fpay* (z-1/3*L) + W_bridge/2*z**2
        z += dz
        v = 2*Fpay + W_bridge * z
        dlist.append(deflection)
        Mlist.append(M)
        vlist.append(v)

    else:
        deflection = -1 / E / Ixx * (Fpay / 6 * (z ** 3 + (z- 1 / 3 * L)**3 + (z- 1 / 2 * L)** 3) + W_bridge / 24 * z ** 4 +A*z + B)
        M = Fpay * z + Fpay * (z - 1 / 3 * L) + Fpay * (z - 1 / 2 * L)+ W_bridge / 2 * z ** 2
        z += dz
        v = 3 * Fpay + W_bridge * z
        dlist.append(deflection)
        Mlist.append(M)
        vlist.append(v)

t= np.linspace(0,L,51)
Mlist=np.array(Mlist)


#Bend stress
Bend_stress= Mlist *0.5*h_bridge/Ixx


p1=plt.figure()
plt.subplot(221)
plt.plot(t,dlist)
plt.xlabel("Length of bridge ")
plt.ylabel("Deflection ")

plt.subplot(222)
plt.plot(t,Mlist)
plt.xlabel("Length of bridge ")
plt.ylabel("Moment ")


plt.subplot(223)
plt.plot(t,Bend_stress)
plt.xlabel("Length of bridge ")
plt.ylabel("Bending stress ")

plt.subplot(224)
plt.plot(t,vlist)
plt.xlabel("Length of bridge ")
plt.ylabel("Shear force ")

plt.show()

#Shear stress
Fshear= 3*Fpay +W_bridge*L
s1= np.linspace(0,w_bridge/2,60)
s2= np.linspace(0,h_bridge,120)
s3= np.linspace(0,w_bridge/2,60)
q0= 1/w_bridge*h_bridge *(-Fshear/Ixx *(thicc_top *h_bridge**2 * w_bridge**2/16 + thicc_side*w_bridge*h_bridge**3/48))
q12= -Fshear/Ixx * thicc_top*h_bridge*s1/2 +q0
q23= -Fshear/Ixx *thicc_side*(h_bridge/2*s2- s2**2/2) -Fshear/Ixx *thicc_top*h_bridge*w_bridge/4 +q0
q34= Fshear/Ixx* thicc_top * h_bridge *s3/2  -Fshear/Ixx * thicc_top*h_bridge*w_bridge/4 +q0

tao12= q12/thicc_top
tao23= q23/thicc_side
tao34= q34/thicc_top

maxtao12= max(tao12,key=abs)
maxtao23= max(tao23,key=abs)
maxtao34= max(tao34,key=abs)

#print(maxtao12, np.where(tao12==maxtao12))
#print(maxtao23, np.where(tao23==maxtao23))
#print(maxtao34, np.where(tao34==maxtao34))

print("Highest shear stress encountered=", max([maxtao12,maxtao23,maxtao34],key=abs)/1E6,'[MPa]')


print('Weight of the bridge=',Volume*rho_bridge,'[kg]')

fig2= plt.figure()
plt.subplot(221)
plt.plot(s1,tao12/1E6)
plt.xlabel("Yes 12")
plt.ylabel("Shearstress [MPa]")

plt.subplot(222)
plt.plot(s2,tao23/1E6)
plt.xlabel("Yes 23")
plt.ylabel("Shearstress [MPa]")

plt.subplot(223)
plt.plot(s3,tao34/1E6)
plt.xlabel("Yes 34")
plt.ylabel("Shearstress [MPa]")

plt.show()



