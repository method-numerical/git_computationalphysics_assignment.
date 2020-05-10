import numpy as np
import matplotlib.pyplot as plt

n=128
fx=np.ones(n)
dx=1
x=np.arange(n)-(n-1)*dx/2

#using numpy to calculate the fourier transform.
fk_np=np.fft.fft(fx)/np.sqrt(n)
k=2*np.pi*np.fft.fftfreq(n) 
ft=dx*np.sqrt(n/2/np.pi)*np.exp(-1j*k*x[0])*fk_np

plt.plot(k,np.real(ft),label='fourier transform of a constant function.')
plt.xlabel('k.')
plt.ylabel('fourier transfrom . .')
plt.legend()
plt.title('fourier transform using the numpy fft.')
plt.show() 