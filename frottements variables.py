import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

duree = 2700
N = 400
g = 9.8


t = np.linspace(0,duree,N) #tableau du temps


def test (Y,y) :

    M = 28.956 #g/mol
    R = 8.314e3
    return np.array([- Y[0] * M * g / (R * Temp(y))])




def Temp(y):
    def aux (y):
        if y < 11e3 : return -6.5*y*1e-3 + 15
        elif y < 20e3 : return -56.5
        elif y < 32e3 : return 1*(y - 20e3)*1e-3 - 56.5
        elif y < 47e3 : return 2.8*(y - 32e3)*1e-3 - 44.5
        elif y < 51e3 : return -2.5
        elif y < 71e3 : return -2.8*(y - 51e3)*1e-3 - 2.5
        else : return -2*(y - 71e3)*1e-3 - 58.5
    return aux(y) + 273.15

def rho_air(y,P):
    y = int(y)
    M_air = 28.956 #! Check this data
    R = 8.314
    return P * M_air / (R * Temp(y))

def k(y,P):
    r1 = 6.4 / 2
    r2 = 8.6 / 2
    S = np.pi * (r1**2 + r2**2)
    Cx = 0.02
    return 1/2 * rho_air(y,P) * Cx * S


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
    #Y[4] => y

    #dv => acceleration
    #dm => dérivée de la masse
    #dy => vitesse selon y
    #dx => vitesse selon x


    if Y[2] > 307240 : #Tant qu'on a du carburant

        dvx = (-F*np.sin(alpha) - k*np.sqrt(Y[0]**2*Y[1]**2)*Y[0])/Y[2]
        dvy = (F*np.cos(alpha) - k*np.sqrt(Y[0]**2*Y[1]**2)*Y[1])/Y[2] - g
        dm = -D
        dy = Y[1]
        dx = Y[0]

    else : #Il n'y a plus de carburant

        dvx = -k*np.sqrt(Y[0]**2*Y[1]**2)*Y[0]/Y[2]
        dvy = - k*np.sqrt(Y[0]**2*Y[1]**2)*Y[1]/Y[2] - g
        dm = 0
        dy = Y[1]
        dx = Y[0]


    if Y[4] < 0 : #On touche le sol

        dvx = 0
        dvy = 0
        dm = 0
        dy = 0
        dx = 0


    return np.array([dvx, dvy, dm, dx, dy])


def crop (t,v) : #On redimensionne les tableaux de temps et de vitesse
    for i in range (len(t)) :
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

# y1 = np.linspace(0,max(y),N)
# P = odeint(test,np.array([1e5]),y1)


#Tracé
# x, yc = crop(x, y)
t, v = crop(t, v)


plt.figure()
# plt.plot(y,P)
plt.xlabel('Altitude en m')
plt.ylabel('Pression en Pa')
plt.figure()
plt.plot(x,y)
plt.xlabel('longitude en m')
plt.ylabel('Altitude en m')

plt.show()