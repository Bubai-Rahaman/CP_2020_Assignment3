'''
Furier transform of a constant function. In these case it is 1.
'''
import numpy as np
import matplotlib.pyplot as plt

numpoints = 1024
xmin = -500
xmax = 500 #xmin and xmax

xarr = np.linspace(xmin,xmax,numpoints)
sampled_data = np.ones(numpoints) #sampling the data
dx = (xmax-xmin)/(numpoints-1)

#calculating dft of the sampled_data
nft = np.fft.fft(sampled_data, norm = 'ortho')
karr = np.fft.fftfreq(numpoints, d=dx) #k-value
karr = 2*np.pi*karr
factor = np.exp(-1j*karr*xmin)

aft = dx*np.sqrt(numpoints/(2*np.pi))*factor*nft



plt.plot(karr,aft,'mo',label = 'Fourier Transform of 1 using np.fft.fft')
plt.xlabel('k')
plt.ylabel('f(k)')
plt.xlim(-.1,.1)
plt.legend()
plt.show()


