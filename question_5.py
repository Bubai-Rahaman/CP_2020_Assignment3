'''
Time comparision to calculated dft using np.fft.fft and my code
'''
import numpy as np
import time
import matplotlib.pyplot as plt

#minimum and maximum value of the data and number of data set
nmin = 4
nmax = 100
nset = 16

dn = int((nmax-nmin)/nset)
narr = np.zeros(nset+1)
my_code_time = np.zeros(nset+1) #array to store time taken by my code
fft_code_time = np.zeros(nset+1) #array to store time taken by np.fft.fft
i = 0
for n in range(nmin,nmax+1,6):
	
	#my code
	start_time_1 = time.time()
	wq = np.zeros(n,dtype = np.complex)
	twq = np.zeros(n, dtype = np.complex)
	T = np.zeros((n,n), dtype = np.complex)
	for q in range(n):
		wq[q] = q
		for p in range(n):
			T[q,p] = np.exp(-1j*2*(np.pi)*q*p/n)/np.sqrt(n)
		twq = T@wq
	end_time_1 = time.time()
	
	#numpy function
	start_time_2 = time.time()
	nft = np.fft.fft(wq, norm = 'ortho')
	end_time_2 = time.time()	
	
	narr[i] = n
	my_code_time[i] = (end_time_1-start_time_1)
	fft_code_time[i] = (end_time_2-start_time_2)
	i = i+1
	
plt.plot(narr,my_code_time, 'b.-',label = 'My code evaluation time')
plt.plot(narr,1000*fft_code_time, 'g.-' ,label = '1000 times of Numpy function evaluation time')
plt.legend()
plt.ylabel('t(sec)')
plt.xlabel('Number of points')
plt.show()	 

