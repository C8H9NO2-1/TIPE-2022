import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

##Champ d'attraction terrestre en fonction de l'altitude
def g(z):
G=6.67e-11
MT=6e24
return G*MT/(z+6400e3)**2

p0=1e5 #Pa
z=np.linspace(0,150000,150000) #altitude

##Température en fonction de l'altitude
def calc_T (z):
if z < 11e3 : return -6.5*z*1e-3 + 15
elif z < 20e3 : return -56.5
elif z < 32e3 : return 1*(z - 20e3)*1e-3 - 56.5
elif z < 47e3 : return 2.8*(z - 32e3)*1e-3 - 44.5
elif z < 51e3 : return -2.5
elif z < 71e3 : return -2.8*(z - 51e3)*1e-3 - 2.5
else : return -2*(z - 71e3)*1e-3 - 58.5

##Équation differentielle de la pression en fonction de l'altitude
def P(Y,z):
M = 28.956 #g/mol
R = 8.314e3
return np.array([-Y[0] * M * g(z) / (R * (calc_T(z) + 273.15))])

##Résolution de l'équation différentielle de la pression
CI = np.array([p0])
p = odeint(P, CI, z)

plt.figure()
plt.plot(z,p, color='red')
plt.ylabel('Pression en Pa')
plt.xlabel('Altitude en m')
plt.grid()
plt.show()
