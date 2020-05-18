import numpy as np
import matplotlib.pyplot as plt

#define sinc function
def f(x):
	try:
		return np.sin(x)/x
	except ZeroDivisionError:
		return 1

#define box function
def tf(k):
	if -1<k<1:
		g = np.sqrt(np.pi/2)
	else:
		g = 0
	return g

#range of sample points
xmin = -500.0
xmax = 500.0 

#importing data from question 2
k2 = np.loadtxt('q_2_fftdata.txt', usecols = 0)
w2 = np.loadtxt('q_2_fftdata.txt', usecols = 1)

#importing data from question 3
k3 = np.loadtxt('q_3_fftdata.txt', usecols = 0)
w3 = np.loadtxt('q_3_fftdata.txt', usecols = 1)


numpoints = 1024
dx = (xmax-xmin)/(numpoints-1)

sampled_data = np.zeros(numpoints)
xarr = np.zeros(numpoints)

#sampling the data
for i in range(numpoints):
	xarr[i] = xmin+i*dx
	sampled_data[i] = f(xarr[i])

#dft of sampled data
nft = np.fft.fft(sampled_data, norm = 'ortho')
karr = np.fft.fftfreq(numpoints, d = dx) 
karr = 2*np.pi*karr	#k-values
factor = np.exp(-1j*karr*xmin) #xmin is non zeros


aft = dx*np.sqrt(numpoints/(2*np.pi))*factor*nft #fourier transform of sinc function
kmin = -10
kmax = 10
dk = (kmax-kmin)/(numpoints-1)
k_plot = np.zeros(numpoints)
tf_plot = np.zeros(numpoints)

#data to plot analytic solution
for i in range(numpoints):
	k_plot[i] = kmin+i*dk
	tf_plot[i] = tf(k_plot[i])

fig = plt.figure()

ax1=fig.add_subplot(121)
ax2 = fig.add_subplot(122)

#ploting the sinc function
ax1.plot(xarr,sampled_data,'g.-', label = 'Sinc function')
ax1.legend()
ax1.set_xlim(-100,100)
ax1.set_xlabel('$x$')
ax1.set_ylabel('$f(x)$')

#ploting the fourier transform of sinc function
ax2.plot(karr,aft,'k*', label = 'FT of sinc function using np.fft.fft')
ax2.plot(k_plot,tf_plot,'m',label = 'Analytically calculated FT of sinc function')
ax2.plot(k2,w2,'g.',label = 'FT of sinc function using FFTW')
ax2.plot(k3,w3,'r.',label = 'FT of sinc function using gsl')
ax2.set_xlabel('k')
ax2.set_ylabel('f(k)')
ax2.legend()

plt.show()
