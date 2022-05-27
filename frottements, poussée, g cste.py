import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import odeint

duree = 2700
N = 400

t = np.linspace(0,duree,N) #tableau du temps


def calc (Y,t) :

    D = 518*9 #Débit massique
    alpha = -10*np.pi/180
    F = 2000e3*9 #Poussée des moteurs
    k = 0.1 #Coefficient de frottements
    g = 9.8

    #Y[0] => Vitesse selon x
    #Y[1] => Vitesse selon y
    #Y[2] => masse m
    #Y[3] => x
    #Y[4] => z

    #dv => acceleration
    #dm => dérivée de la masse
    #dz => vitesse selon y
    #dx => vitesse selon x


    if Y[2] > 307240 : #Tant qu'on a du carburant 

        dvx = (-F*np.sin(alpha) - k*np.sqrt(Y[0]**2*Y[1]**2)*Y[0])/Y[2] 
        dvy = (F*np.cos(alpha) - k*np.sqrt(Y[0]**2*Y[1]**2)*Y[1])/Y[2] - g
        dm = -D
        dz = Y[1]
        dx = Y[0]

    else : #Il n'y a plus de carburant 

        dvx = -k*np.sqrt(Y[0]**2*Y[1]**2)*Y[0]/Y[2]
        dvy = - k*np.sqrt(Y[0]**2*Y[1]**2)*Y[1]/Y[2] - g
        dm = 0
        dz = Y[1]
        dx = Y[0]

    
    if Y[4] < 0 : #On touche le sol

        dvx = 0
        dvy = 0
        dm = 0
        dz = 0
        dx = 0


    return np.array([dvx, dvy, dm, dx, dz])


def crop (t,v) : #On redimensionne les tableaux de temps et de vitesse
    for i in range (len(t)+1) : 
        if v[i]==v[i+1] : return t[:i], v[:i]


Y0 = np.array([0,0,1807240,0,0]) #Conditions initiales
sol = odeint(calc,Y0,t) #résolution de(s) équa. diff.
n = len(sol)

#Extraction des données
vx = np.array([sol[i][0] for i in range (n)])
vy = np.array([sol[i][1] for i in range (n)])
v = np.sqrt(vx**2+vy**2)
x = np.array ([sol[i][3] for i in range (n)])
y = np.array ([sol[i][4] for i in range (n)])
m = np.array ([sol[i][2] for i in range (n)])

#Tracé
plt.figure()
plt.plot(x,y)
plt.figure()
plt.plot(crop(t,v)[0],crop(t,v)[1])
plt.show()

