import numpy as np
import matplotlib.pyplot as plt

#Analytic Fourier transform of Gaussian
def fun(k):
	return np.exp(-k**2/4)/np.sqrt(2)

#importing data
k = np.loadtxt('q_4_fftdata.txt', usecols = 0)
w = np.loadtxt('q_4_fftdata.txt', usecols = 1)

k_plot = np.linspace(-40,40,201)
tf = fun(k_plot)

plt.plot(k_plot,tf,'g',label = "Analytically calulated FT of Gaussian")
plt.plot(k,w,'r.', label = 'Numerically calculated(FFTW) FT of Gaussian')
plt.legend()
plt.xlabel('k')
plt.ylabel('tf')
plt.xlim(-20,20)
plt.show()
