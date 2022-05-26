import numpy as np
import matplotlib.pyplot as plt 

def masse (t,nb_moteurs=9): #masse en fonction du temps

    n = 518*t*nb_moteurs
    
    if n > 1500000:
        n = 1500000
        return 1807240 - n, True
    else:
        return 1807240 - n, False #kg

def propulsion (out_of_fuel,nb_moteurs=9):
    if out_of_fuel:
        return 0
    else:
        return 1961e3 * nb_moteurs #N


def trajectoire1 (duree, angle):
    g = 9.8

    alpha = -angle*np.pi/180
    
    N = 400
    
    x0 = 0
    y0 = 0
    v0 = 0
    
    m, out_of_fuel = masse(0)
    F = propulsion(out_of_fuel)


    X = np.zeros(N)
    Y = np.zeros(N)
    V = np.zeros(N)
    
    # Initialisation des tableaux
    X[0] = x0
    Y[0] = y0
    V[0] = v0
    
    T = np.linspace(0,duree,N) #tableau numpy du temps

    for i in range (1,N) : 

        if out_of_fuel : #Calculs lorsque le carburant est consommé

            t_1 = T[i] - t_0 #Nouvelle origine du temps

            #Calcul de la vitesse 
            vx_1 = vx_0 
            vy_1 = -g*t_1 + vy_0

            #Calcul de la position
            x_1 = 1/2 * vx_1 * t_1 + x_0
            y_1 = 1/2 * vy_1 * t_1 + y_0

            #Enregistrement dans le tableau
            X[i] = x_1
            Y[i] = y_1
            V[i] = np.sqrt(vx_1**2 + vy_1**2)

        else : #Calculs lorque le carburant est n'est pas encore consommé

            t_0 = T[i]

            #Calcul de la vitesse
            vx_0 = -F/m*np.sin(alpha)*t_0
            vy_0 = (F/m*np.cos(alpha)-g)*t_0

            #Calcul de la position
            x_0 = 1/2 * vx_0 * t_0
            y_0 = 1/2 * vy_0 * t_0

            #Enregistrement dans le tableau
            X[i] = x_0
            Y[i] = y_0
            V[i] = np.sqrt(vx_0**2 + vy_0**2)

            m, out_of_fuel = masse(T[i]) #Calcul de la masse
            F = propulsion(out_of_fuel) #Calcul de la propulsion

        if Y[i] < 0 and i > 10 : #Si la fusée touche le sol on stop

            return X[:i],Y[:i],V[:i],T[:i]

    return X,Y,V,T

x, y, v, t = trajectoire1(2700,10)

#Graphe de la trajectoire
plt.figure() 
plt.plot(x,y)
plt.xlabel ('x en m')
plt.ylabel ('z en m')
#Graphe de la vitesse en fonction du temps
plt.figure()
plt.plot(t,v)
plt.xlabel ('t en sec')
plt.ylabel ('vitesse en m/s')

plt.show()