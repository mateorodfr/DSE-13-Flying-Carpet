import numpy as np

h_max= 400
c= 299792458
r= 400/np.cos(1/3*np.pi)
k=1.38*10**(-23)
L_cable1= 1
L_cable2=4
Bitrate= 25*10**6
Bandwith= (54-30)*10**6
D= 1



#Frequency
Freq_max= 54E6
Wavelength_max= c/Freq_max

#Antenna gains
G_receiver= 2.14
G_trans = 10* np.log10(np.pi**2 * D**2/Wavelength_max**2 )


#Noise
L1= 0.5904 * L_cable1
L2= 0.5904 * L_cable2
G_amplifier= 20
G_cable1= 1/L1
G_cable2= 1/L2
F_amp= 10**(5.7/10)
F_receiver= 10**(1.6/10) #db
T0= 290
T_amp= T0 * (F_amp -1)
T_antenna= 1005000
T_cable1= T0 *((1-L1)/L1)
T_cable2= T0 *((1-L2)/L2)
T_receiver= T0 * (F_receiver-1)
Tsys= T_antenna + T_cable1/G_trans + T_amp/G_trans/G_cable1 + T_cable2/G_amplifier/G_trans/G_cable1 +T_receiver/ G_amplifier/ G_trans/ G_cable2/G_cable1

Noise0= k*Tsys


#Loss factors
alpha= 5
VSWR=1.5
Lfs= -10* alpha *np.log10(4*np.pi*r/Wavelength_max) #Urban env
#Lfs= 10 * np.log10((Wavelength_max/(4*np.pi*r))**2) #Normal env
#La=  -0.5904 * L_cable1 -0.5904 * L_cable2  #Transmission path loss Read it off the graph
Ll= -0.5904 * L_cable1  #Loss factor going to antenna
Lpr= ... #Pointing error loss
Lr= -14 # Any other losses that may occur

#Power & gains
P_trans_max= 10*np.log10(250)

P_receiver_max= P_trans_max + G_receiver+ G_trans + Lfs + Ll +Lr #+ Lpr
PNR= 10**(P_receiver_max/10)/Noise0
SNR= 10**(P_receiver_max/10)/Noise0*Bitrate
Capacity= Bandwith* np.log2(1+SNR)

#print(Lfs)
#print(La)
print("Noise",Noise0)
print(P_receiver_max)
print("Power Received", 10**(P_receiver_max/10))

print("PNR=", 10*np.log10(PNR))
print("SNR ratio=", SNR)
print("SNR db=", 10*np.log10(SNR))
print("Capacity=", Capacity/10**6, '[Mbit/s]')

print(G_trans)
print(Wavelength_max)


