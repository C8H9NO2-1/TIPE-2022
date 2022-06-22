import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

duree = 3000
N = 800
g = 9.8


t = np.linspace(0,duree,N) #tableau du temps

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

def angle_phase1(y):
    def aux(y):
        if y < 10e3: return -5-5*y/10e3
        elif y < 30e3: return -10 - 20 * (y - 10e3) / 20e3
        elif y < 60e3: return -30 - 40 * (y - 30e3) / 30e3
        else: return -70 - 20 * (y - 60e3) / 20e3
    return aux(y) * np.pi / 180

# aux(y) * np.pi / 180
def calc (Y,t) :

    D = 518*9 #Débit massique
    alpha = angle_phase1(Y[4])
    # print(alpha * 180 / np.pi) 
    F = 2000e3*9 #Poussée des moteurs
    
    g = 9.8
    M = 28.956e-3 #g/mol
    R = 8.314

    #Coefficient de frottements
    r1 = 6.4 / 2
    r2 = 8.6 / 2
    S = np.pi * (r1**2 + r2**2)
    Cx = 0.02
    k = Y[5] * M / (R * Temp(Y[4])) * Cx * S


    #?Y[0] => Vitesse selon x
    #?Y[1] => Vitesse selon y
    #?Y[2] => masse m
    #?Y[3] => x
    #?Y[4] => y
    #?Y[5] => Pression P

    #?dv => acceleration
    #?dm => dérivée de la masse
    #?dy => vitesse selon y
    #?dx => vitesse selon x
    #?dP => dérivée de la Pression

    if Y[2] > 307240 : #Tant qu'on a du carburant

        dvx = (-F*np.sin(alpha) - k*np.sqrt(Y[0]**2 + Y[1]**2)*Y[0])/Y[2]
        dvy = (F*np.cos(alpha) - k*np.sqrt(Y[0]**2 + Y[1]**2)*Y[1])/Y[2] - g
        dm = -D
        dy = Y[1]
        dx = Y[0]
        dP = (-Y[5] * M * g / (R * Temp(Y[4]))) * abs(Y[1])
        

    else : #Il n'y a plus de carburant
    
        dvx = -k*np.sqrt(Y[0]**2 + Y[1]**2)*Y[0]/Y[2]
        dvy = - k*np.sqrt(Y[0]**2 + Y[1]**2)*Y[1]/Y[2] - g
        dm = 0
        dy = Y[1]
        dx = Y[0]
        dP = (-Y[5] * M * g / (R * Temp(Y[4]))) * abs(Y[1])

    if Y[4] < 0 : #On touche le sol

        dvx = 0
        dvy = 0
        dm = 0
        dy = 0
        dx = 0
        dP = 0
        

    return np.array([dvx, dvy, dm, dx, dy, dP])


def crop (t,v,x,y,P) : #On redimensionne les tableaux lorsque les valeurs deviennent constantes
    for i in range (len(t)-1) :
        if y[i]==y[i+1] : return t[:i], v[:i], x[:i], y[:i], P[:i]


Y0 = np.array([0,0,1807240,0,0,1e5]) #Conditions initiales
sol = odeint(calc,Y0,t) #résolution de(s) équa. diff.
n = len(sol)

#Extraction des données
vx = np.array([sol[i][0] for i in range (n)])
vy = np.array([sol[i][1] for i in range (n)])
v = np.sqrt(vx**2+vy**2)
x = np.array ([sol[i][3] for i in range (n)])
y = np.array ([sol[i][4] for i in range (n)])
m = np.array ([sol[i][2] for i in range (n)])
P = np.array ([sol[i][5] for i in range (n)])

t,v,x,y,P = crop(t,v,x,y,P)


#Tracé

plt.figure()
plt.plot(x[:-1], y[:-1], marker='x', ls='none')
plt.ylabel('Altitude en m')
plt.xlabel('Longueur en m')
plt.grid()
plt.show()
