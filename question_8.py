'''
2d Fourier transform of gaussian function
'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#define Gaussian fuction
def Gaussian_2d(x,y):
  return np.exp(-(x**2+y**2)) 	

#Analytical FT of Gaussian  
def fft_Gaussian_2d(kx,ky):
  return np.exp(-(kx**2+ky**2)/4)/2

#setting x and y range
Numpoints = 100
xmin,xmax = -20, 20
ymin,ymax = -20, 20

dx=(xmax-xmin)/(Numpoints-1)
dy=(ymax-ymin)/(Numpoints-1)

xarr = np.linspace(xmin,xmax,Numpoints, dtype = np.complex)
yarr = np.linspace(ymin,ymax,Numpoints, dtype = np.complex)

sampled_data = np.zeros([Numpoints,Numpoints], dtype=np.complex_)

#kx and ky values
kx = 2*np.pi*np.fft.fftfreq(Numpoints,dx)
ky = 2*np.pi*np.fft.fftfreq(Numpoints,dy)

#sampling the data
for i in range(Numpoints):
  for j in range(Numpoints):
    sampled_data[i][j] = Gaussian_2d(xarr[i],yarr[j])
    
    
nft = np.fft.fft2(sampled_data, norm='ortho') #dft of the sampled data

aft = np.zeros([Numpoints,Numpoints], dtype=np.complex_)

for i in range(Numpoints):
  for j in range(Numpoints):
    aft[i][j]=dx*dy*(Numpoints/(2.0*np.pi))*(np.exp(-1j*kx[i]*xmin - 1j*ky[j]*ymin ))*nft[i][j]

fig = plt.figure()
ax = plt.axes(projection="3d")

x = np.zeros(Numpoints*Numpoints)
y = np.zeros(Numpoints*Numpoints)
Fk1 = np.zeros(Numpoints*Numpoints)
Fk2 = np.zeros(Numpoints*Numpoints)
for i in range(Numpoints):
  for j in range(Numpoints):
    x[j+i*Numpoints] = kx[i]
    y[j+i*Numpoints] = ky[j]
    Fk1[j+i*Numpoints] = aft[i][j].real #storing Numerical result
    Fk2[j+i*Numpoints] = fft_Gaussian_2d(kx[i],ky[j]) #storing analytical result
   
ax.plot3D(x, y, Fk2, 'b.',label="Analytic Fourier Transform of Gaussian")
ax.plot3D(x, y, Fk1, 'r.',label="Numerical Fourier Transform of Gaussian(np.fft.fft2)")
ax.set_xlabel('kx')
ax.set_ylabel('ky')
ax.set_zlabel('F(kx,ky)')
plt.legend()
plt.show()

