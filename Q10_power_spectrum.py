#connect to internet, because this code requireinternet access.

import numpy as np
import matplotlib.pyplot as plt

f1='http://theory.tifr.res.in/~kulkarni/noise.txt'
data=np.loadtxt(f1,usecols=0)

n=len(data)
x=np.arange(n)

dft=np.fft.fft(data)
k=2*np.pi/n*(np.arange(n)-n/2) 
l=np.argsort(k)
dft=dft[l]

periodogram=abs(dft)*abs(dft)/n

bin_length=n//10
extra_elements=n%10
bin_edge=[]
pow_spec_av=[]
for i in range(0,9):
  bin_index=np.arange(bin_length)+i*bin_length
  bin_edge.append(k[bin_index[0]])
  pow_spec_av.append(np.mean(periodogram[bin_index]))

last_bin_index=np.arange(9*bin_length,n)
bin_edge.append(k[9*bin_length])
pow_spec_av.append(np.mean(periodogram[last_bin_index]))

plt.figure(figsize=[7,20])

plt.subplot(411)
plt.plot(x,data,label='measurement.')
plt.xlabel('x.')
plt.legend()

plt.subplot(412)
plt.plot(k,np.real(dft),label='dft.')
plt.xlabel('k.')
plt.legend()

plt.subplot(413)
plt.plot(k,periodogram,label='power spectrum using periodogram....')
plt.xlabel('k.')
plt.legend()

plt.subplot(414)
plt.bar(bin_edge,pow_spec_av,align='edge',width=1,label='binned power spectrum .')
plt.xlabel('k.')
plt.legend()

plt.show()


