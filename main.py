import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

##Champ d'attraction terrestre en fonction de l'altitude
def g(z):
	G=6.67e-11
	MT=6e24
	return G*MT/(z+6400e3)**2

p0=1e5 #Pa
z=np.linspace(0,1500000,1500000) #altitude

##Température en fonction de l'altitude
def Temp (z):
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
	return np.array([-Y[0] * M * g(z) / (R * (Temp(z) + 273.15))])

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
    n = 518*t*nb_moteurs
    
    if n > 1500000:
        return 1807240, True
    else:
        return 1807240 - n, False #kg
    
def propulsion (out_of_fuel,nb_moteurs=9):
    if out_of_fuel:
        return 0
    else:
        return 1961e3 * nb_moteurs #N

def rho_air(z):
    y = round(z)
    M_air = 0.01801528 #! Check this data
    R = 8.314
    return p[y][0] * M_air / (R * (Temp(z) + 273.15))

def frottements(z, m, v):
    r1 = 6.4 / 2
    r2 = 8.6 / 2
    S = np.pi * (r1**2 + r2**2)
    Cx = 0.02
    return 1/2 * rho_air(z) * Cx * S * v**2


def trajectoire ():
    N = 400
    
    x0 = 0
    y0 = 0
    
    out_of_fuel = False
    
    X = np.zeros(N)
    Y = np.zeros(N)
    
    X[0] = x0
    Y[0] = y0
    
    v = 0
    
    T = np.linspace(0,60,N) #tableau numpy du temps
    
    for i in range (1,N) : 
    
        m, out_of_fuel = masse(T[i])
        p = propulsion(out_of_fuel)
        # X[i] = 1/2*g(Y[i-1])*np.sin(alpha)*T[i]**2
        v = (p/m - g(Y[i-1]) - frottements(Y[i - 1], m, v) * np.sign(v)) * T[i]
        Y[i] = 1/2*v*T[i]
        
        
    return X,Y
    

plt.plot(trajectoire()[0],trajectoire()[1])
plt.xlabel ('x en m')
plt.ylabel ('z en m')
plt.show()
