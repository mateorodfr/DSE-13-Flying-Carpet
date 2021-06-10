import numpy as np

h_max= 400
c= 299792458
r= 400/np.cos(1/3*np.pi)
k=1.38E-23
L_cable= 5
Bitrate= 1/128E3
Bandwith= (54-30)*10**6

G_receiver= 2.14
G_trans = 2.14

#Frequency
Freq_max= 54E6
Wavelength_max= c/Freq_max


#Noise
L= 0.5904 * L_cable
G_amplifier= 30
G_cable= 1/L
F_amp= 10**(4/10)
F_receiver= 10**(1.6/10) #db
T0= 290
T_amp= T0 * (F_amp -1)
T_antenna= 1005000
T_cable= T0 *((1-L)/L)
T_receiver= T0 * (F_receiver-1)
Tsys= T_antenna + T_amp/G_receiver + T_cable/G_amplifier/G_trans +T_receiver/ G_amplifier/ G_trans/ G_cable


Noise0= k*Tsys



#Loss factors
alpha= 5
Lfs= -10* alpha *np.log10(4*np.pi*r/Wavelength_max) #Urban env
#Lfs= 10 * np.log10((Wavelength_max/(4*np.pi*r))**2) #Normal env
La=  0.5904 * L_cable  #Transmission path loss
#Ll= ...  #Loss factor antenna
#Lpr= ... #Pointing error loss
#Lr= ...  # Any other losses that may occur

#Power & gains
P_trans_max= 10*np.log10(250)



P_receiver_max= P_trans_max + G_receiver+ G_trans + Lfs + La
PNR= 10**(P_receiver_max/10)/Noise0
SNR= 10**(P_receiver_max/10)/Noise0/Bitrate
Capacity= Bandwith* np.log2(1+SNR)

#print(Lfs)
#print(La)
print("Noise",Noise0)
#print(P_receiver_max)
print("Power Received", 10**(P_receiver_max/10))

print("PNR=", PNR)
print("SNR=", SNR)
print("Capacity=", Capacity/10**6, '[Mbit/s]')

