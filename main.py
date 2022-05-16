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

#plt.figure()
#plt.plot(z,p, color='red')
#plt.ylabel('Pression en Pa')
#plt.xlabel('Altitude en m')
#plt.grid()
#plt.show()

##Première phase du vol

nb_moteurs = 9 #nombre de moteurs en marche

alpha = 179.5*np.pi/180 #rad

def masse (t,nb_moteurs=9): #masse en fonction du temps
    return 1807240 - 518*t*nb_moteurs #kg
    
def propulsion (nb_moteurs=9):
    
    return 1961e3 * nb_moteurs #N
    


def trajectoire ():
    N = 400
    
    x0 = 0
    y0 = 0
    
    X = np.zeros(N)
    Y = np.zeros(N)
    
    X[0] = x0
    Y[0] = y0
    
    T = np.linspace(0,400,N) #tableau numpy du temps
    
    for i in range (1,N) : 
    
        if masse(T[i]) < 0 : 
            X[i] = 1/2*g(Y[i-1])*np.sin(alpha)*T[i]**2
            Y[i] = 1/2*(propulsion()/masse(T[i])-g(Y[i-1])*np.cos(alpha))*T[i]**2
        
        
    return X,Y
    

plt.plot(trajectoire()[0],trajectoire()[1])
plt.xlabel ('x en m')
plt.ylabel ('z en m')
plt.show()
