import numpy as np

h_max= 400
c= 3 * 10**8
r= 400/np.cos(1/3*np.pi)
k=1.38*10**(-23)
L_cable1= 1
L_cable2=4
Bitrate= 12*10**6
Bandwith= (54-30)*10**6
D = 1.5



#Frequency
Freq_max= 54*10**6     # We need a higher frequency antenna, else the receiving end is going to be weak, making it lossy.
# However, going higher in frequency increases the L_FS, thus making it worse
Wavelength_min= c/Freq_max
# G_high = 1
# A_required = G_high / 0.9 * (Wavelength_min**2/4/np.pi)
# print(A_required)
#Antenna gains
G_trans_dB= 2.14
G_trans = 10**(G_trans_dB/10)
G_receiver_dB = 10* np.log10(np.pi**2 * D**2/Wavelength_min**2)
Trans_eff = 0.9
G_receiver = np.pi * D**2/Wavelength_min**2 * Trans_eff


#Noise
L1_dB= 0.5904 * L_cable1
L1 = 10**(L1_dB/10)
L2_dB= 0.5904 * L_cable2
L2 = 10**(L2_dB/10)
G_amplifier_dB= 20
G_amplifier = 10**(G_amplifier_dB/10)
G_cable1= 1/L1
G_cable2= 1/L2
F_amp= 10**(5.7/10)
F_receiver= 10**(1.6/10) # NOT dB
T0= 290
T_amp= T0 * (F_amp-1)
T_antenna= 750
T_cable1= T0 *(L1-1)
T_cable2= T0 *(L2-1)
T_receiver= T0 * (F_receiver-1)
Tsys= T_antenna + T_cable1 + T_amp/G_cable1 + T_cable2/G_amplifier/G_cable1 +T_receiver/ G_amplifier/ G_cable2/G_cable1


B_noise = Bandwith * 6.28/4  # The actual value needs to be computed by creating a filter!
Noise0= k*Tsys
N_tot = Noise0 * B_noise

#Loss factors
alpha= 3.6    # Okey, this Alpha is waaaaay too high. We need either a lot stronger amplifier or larger antenna
#VSWR=1.5
Lfs_dB= -10* alpha *np.log10(4*np.pi*r/Wavelength_min) #Urban env
Lfs = 10**(Lfs_dB/10)
#Lfs= 10 * np.log10((Wavelength_max/(4*np.pi*r))**2) #Normal env
#La=  -0.5904 * L_cable1 -0.5904 * L_cable2  #Transmission path loss Read it off the graph
Ll_dB= -0.5904 * L_cable1  #Loss factor going to antenna
Ll = 10 ** (Ll_dB/10)
Lpr= ... #Pointing error loss
#Lr_dB= -14 # Any other losses that may occur
#Lr = 10**(Lr_dB/10)


#Power & gains
P_trans_max = 250
P_trans_max_db= 10*np.log10(P_trans_max)

P_receiver_max= P_trans_max * G_receiver * G_trans * Lfs * Ll #* Lr #* Lpr
#PNR= P_receiver_max/N_tot
SNR= P_receiver_max/N_tot
SNR_1 = P_receiver_max/Noise0/Bitrate
Capacity= Bandwith * np.log2(1+SNR_1)

#print(Lfs)
#print(La)
print("Noise",10*np.log10(Noise0))
print("Power Received", P_receiver_max)

#print("PNR=", PNR)
print("SNR_1 ratio=", 10*np.log10(SNR_1))
print("SNR =", 10*np.log10(SNR))
print("Capacity=", Capacity/10**6, '[Mbit/s]')

print(G_receiver)



