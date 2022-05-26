import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

##Champ d'attraction terrestre en fonction de l'altitude
def g(z):
	G=6.67e-11
	MT=6e24
	return G*MT/(z+6400e3)**2

p0=1e5 #Pa
z=np.linspace(0,15000000,15000000) #altitude

##Température en fonction de l'altitude en Kelvin
def Temp(z):
    def aux (z):
        if z < 11e3 : return -6.5*z*1e-3 + 15
        elif z < 20e3 : return -56.5
        elif z < 32e3 : return 1*(z - 20e3)*1e-3 - 56.5
        elif z < 47e3 : return 2.8*(z - 32e3)*1e-3 - 44.5
        elif z < 51e3 : return -2.5
        elif z < 71e3 : return -2.8*(z - 51e3)*1e-3 - 2.5
        else : return -2*(z - 71e3)*1e-3 - 58.5
    return aux(z) + 273.15

##Équation differentielle de la pression en fonction de l'altitude
def P(Y,z):
	M = 28.956 #g/mol
	R = 8.314e3
	return np.array([-Y[0] * M * g(z) / (R * Temp(z))])

##Résolution de l'équation différentielle de la pression
CI = np.array([p0])
p = odeint(P, CI, z)

# plt.figure()
# plt.plot(z,p, color='red')
# plt.ylabel('Pression en Pa')
# plt.xlabel('Altitude en m')
# plt.grid()
# plt.show()

##Première phase du vol

nb_moteurs = 9 #nombre de moteurs en marche

alpha = -50*np.pi/180 #rad

def masse_ensemble(t,nb_moteurs=9): #masse de l'ensemble booster avion en fonction du temps
    n = 518*t*nb_moteurs
    
    if n > 1500000:
        return 328000, True
    else:
        return 1807240 - n, False #kg
    
def propulsion(out_of_fuel,nb_moteurs=9):
    if out_of_fuel:
        return 0
    else:
        return 1961e3 * nb_moteurs #N

def rho_air(z):
    y = round(z)
    M_air = 0.01801528 #! Check this data
    R = 8.314
    return p0 * M_air / (R * Temp(z))

def frottements(z, v):
    r1 = 6.4 / 2
    r2 = 8.6 / 2
    S = np.pi * (r1**2 + r2**2)
    Cx = 0.02
    return 1/2 * rho_air(z) * Cx * S * v**2 * np.sign(-v)


def trajectoire1 ():
    N = 400
    
    x0 = y0 = 0
    out_of_fuel = False
    
    X = np.zeros(N)
    Y = np.zeros(N)
    X[0] = x0
    Y[0] = y0
    
    vx_1 = vy_1 = v_1 = 0
    
    fx = fy = frottements(0, 0)
    
    T = np.linspace(0,400,N) #tableau numpy du temps
    
    for i in range (1,N) : 
    
        m, out_of_fuel = masse_ensemble(T[i])
        p = propulsion(out_of_fuel)
        
        #* Selon la base 1:
        vx_1 = (-np.sin(alpha) * g(Y[i - 1]) + fx / m) * T[i]
        vy_1 = (-np.cos(alpha) * g(Y[i - 1]) + (fy + p) / m) * T[i]
        v_1 = np.sqrt(vx_1**2 + vy_1**2)
        
        #* Selon la base 0:
        vx_0 = np.cos(alpha) * vx_1 - np.sin(alpha) * vy_1
        vy_0 = np.sin(alpha) * vx_1 + np.cos(alpha) * vy_1
        
        # print(vy_0)
        
        X[i] = 1/2 * vx_0 * T[i]
        Y[i] = 1/2*vy_0*T[i]
        
        if Y[i] >= 80e3:
            print(m - 328000)
            return X, Y, T, i
        
        #? Comme ça la fusée ne s'enfonce pas dans le sol
        if Y[i] < 0:
            Y[i] = 0
            X[i] = X[i - 1]
            # vx_1 = vy_1 = v_1 = vx_0 = vy_0 = 0
            
        fx = frottements(Y[i], vx_1)
        fy = frottements(Y[i], vy_1)
        
        print("Altitude // ", i, ": ", Y[i])
        
    return X,Y,T, N - 1 # This line is just here to make the code prettier
    
def trajectoire2(x0, y0):
    N = 3600
    
    X = np.zeros(N)
    Y = np.zeros(N)
    X[0] = x0
    Y[0] = y0
    
    v = 0 #! Ici aussi
    
    T = np.linspace(0, 3600, N) #tableau numpy du temps
    
    for i in range(1, N):
        pass

    return X, Y

x, y, t, i = trajectoire1()
x = x[:i]
y = y[:i]
plt.plot(x,y)
# plt.xlabel ('x en m')
# plt.ylabel ('z en m')
plt.show()
