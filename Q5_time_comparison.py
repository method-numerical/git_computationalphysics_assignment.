import numpy as np
import matplotlib.pyplot as plt
import time as t 

def dft(f):
  n=len(f)
  dft=[]
  for q in range(n):
    a=0
    for p in range(n):
      a+=f[p]*np.exp(-2j*np.pi*p*q/n)
    dft.append(a)
  return dft/np.sqrt(n)

sample=np.arange(4,100,12)
t_code=[]
t_numpy=[]

for i in sample:
  arr=np.arange(i)
  
  t1=t.perf_counter()
  sol1=dft(arr)
  t2=t.perf_counter()
  t_code.append(10**9*(t2-t1))
  
  t3=t.perf_counter()
  sol2=np.fft.fft(arr)/np.sqrt(i)
  t4=t.perf_counter()
  t_numpy.append(10**9*(t4-t3))

plt.plot(sample,t_code,label='time for the code ~ n².')
plt.plot(sample,t_numpy,label='time for numpy fft ~ n*ln(n)/ln2.')
plt.xlabel('n.')
plt.ylabel('t.')
plt.title('comparison between general dft and fft algorithm. ')
plt.legend()
plt.show()
  
      