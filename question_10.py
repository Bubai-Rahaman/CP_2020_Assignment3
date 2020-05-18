'''
Power Spectrum
'''

import numpy as np
import matplotlib.pyplot as plt

#importing data
noise = np.loadtxt("noise.txt", usecols = 0)
t = np.arange(0,noise.size,1)

#calculating dft of the noise
dft_noise = np.fft.fft(noise, norm = 'ortho')
omega = np.fft.fftfreq(t.size, d =1)

#calculating power spectrum
p =dft_noise*np.conj(dft_noise)/dft_noise.size

bins_data = np.zeros(10)
bin_number = np.arange(1,bins_data.size+1,1)
for i in range(bins_data.size):
	s = 0
	for j in range(int(i*51),int(51*(i+1)),1):
		s = s + p[j]
	bins_data[i] = s/51



#ploting noise
plt.figure()
plt.plot(t,noise,'r.-',label = 'noise_data')
plt.xlabel('t')
plt.ylabel('Noise')
plt.legend()

#ploting dft of the noise
plt.figure()
plt.plot(omega, dft_noise,'m^', label = 'DFT of noise')
plt.xlabel('frequency')
plt.ylabel('DFT of noise')
plt.legend()

#ploting power spectrum
plt.figure()
plt.plot(omega,p,'g.', label = 'Power spectrum of the noise')
plt.xlabel('frequency')
plt.ylabel('Power spectrum')
plt.legend()

#bin power spectrum
plt.figure()
plt.plot(bin_number,bins_data,'b*',label = 'Bin power spectrum')
plt.xlabel('Bin number')
plt.ylabel('Power spectrum')
plt.legend()

plt.show()
