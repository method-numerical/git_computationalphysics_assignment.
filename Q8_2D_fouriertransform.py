import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm 
from mpl_toolkits.mplot3d import Axes3D

n=128
dx=0.7
xmin=-(n-1)*dx/2
xmax=xmin+(n-1)*dx 

x=np.linspace(xmin,xmax,n)
y=x
[X,Y]=np.meshgrid(x,y)

Z=np.exp(-X**2-Y**2) 

kx=2*np.pi/(n*dx)*(np.arange(n)-n/2)
ky=kx
[KX,KY]=np.meshgrid(kx,ky)

dft=np.fft.fft2(Z)/n 

#sorting dft.
k=np.fft.fftfreq(n,dx)
l=np.argsort(k)
for i in range(n):
  dft[i]=dft[i][l]
for i in range(n):
  dft[:,i]=dft[:,i][l]
  
dft=(dx*np.sqrt(n/(2*np.pi)))**2*np.exp(-1j*(KX+KY)*xmin)*dft 

#analytic ft..
Z_an=np.exp(-(KX**2+KY**2)/4)/2

fig=plt.figure()

ax=fig.add_subplot(211,projection='3d')
ax.plot_surface(KX,KY,np.real(dft),cmap=cm.coolwarm)
ax.set_xlabel('kx..')
ax.set_ylabel('ky..')
ax.set_zlabel('ft(kx,ky).')
ax.set_title("numerical fourier transform of 2D gaussian.") 

ax=fig.add_subplot(212,projection='3d')
ax.plot_surface(KX,KY,Z_an,cmap=cm.jet,label='analytic ft..')
ax.set_xlabel('kx..')
ax.set_ylabel('ky..')
ax.set_zlabel('ft(kx,ky).')
ax.set_title("analytical fourier transform of 2D gaussian.")

plt.show()


















 
