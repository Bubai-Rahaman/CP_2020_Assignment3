'''
Convolution of box function with itself
'''

import numpy as np
import matplotlib.pyplot as plt

#define box function
def box(x):
	if -1<x<1:
		return 1
	else:
		return 0

xmin,xmax = -10,10 #x-range
N = 512
dx = (xmax-xmin)/(N-1)

xarr = np.zeros(N)
sampled_data = np.zeros(N)

for i in range(N):
	xarr[i] = xmin + i*dx
	sampled_data[i] = box(xarr[i])

dft = np.fft.fft(sampled_data, norm = 'ortho') #dft of the sampled data
fx = np.fft.ifft(dft*dft, norm = 'ortho') 
fx = dx*np.sqrt(N)*fx #convolution function 

shift_fx = np.fft.fftshift(fx) #shifiting zero frequency part to the centre 

plt.plot(xarr,shift_fx,'b.', label = 'Convolution function of box function with itself')
plt.plot(xarr,sampled_data, 'k', label = 'Box function')
plt.xlabel('x');
plt.ylabel('y(x)')
plt.xlim(-5,5)
plt.legend()
plt.show()
