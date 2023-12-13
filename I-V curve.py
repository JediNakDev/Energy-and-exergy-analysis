from power import findI
import numpy as np
import matplotlib.pyplot as plt

G = [1000,1000,1000,1000,1000]
V = np.arange(-5,55,0.1)
T = [273.15+10,273.15+25,273.15+40,273.15+55,273.15+70]

I = findI(V,T,G)
P=[]
for i in range(5):
    P.append(V * I[i])
plt.plot(V,P[0])
plt.plot(V,P[1])
plt.plot(V,P[2])
plt.plot(V,P[3])
plt.plot(V,P[4])
plt.legend(['10°C','25°C','40°C','55°C','70°C'])
plt.xlabel('Voltage (V)')
plt.ylabel('Power (W)')
plt.xlim(-5,55)
plt.ylim(0,400)
plt.minorticks_on()
plt.show()