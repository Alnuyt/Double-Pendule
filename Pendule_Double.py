# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 14:51:49 2022

@author: alexandre
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
plt.rcParams['animation.writer'] = "avconv"

#definition des paramètres
N = 60000
t_max = 60
t = np.linspace(0, t_max, N)
h = t[1]
g = 9.81
l = 1
m = 1
omega = np.sqrt(g/l)


#implementation de la fonction
#fonction qui prend en argument l'etat du systeme au temps t et renvoie sa derivee
def derivative_1(y):
    theta, phi, thetadot, phidot = y
    s = np.sin(theta - phi)
    c = np.cos(theta - phi)
    dthetadot = ((((-1) * s * (((thetadot)**2 * c) + phidot**2)) + omega**2 * ((np.sin(phi) * c) - 2 * np.sin(theta)))/(1 + s**2))
    dphidot = (( s * (c * (phidot)**2 + 2 * thetadot**2) + 2 * omega**2 * ((np.sin(theta) * c) - np.sin(phi)))/(1 + s**2))
    dydt = np.array([thetadot, phidot, dthetadot, dphidot])
    return dydt

#creation matrice solution + conditions initiales
c_1 = [(np.pi)/ 20, (np.pi)/ 30, 0, 0]
c_2 = [(np.pi)/ 20, (np.pi)/ 30 + (np.pi)/ 60, 0, 0]
c_3 = [(np.pi)/ (2), (np.pi)/2, 0, 0]
c_4 = [(np.pi)/(2), (np.pi)/2 + (np.pi)/30, 0, 0]
sol_1 = np.zeros((4, N))
sol_2 = np.zeros((4, N))
sol_3 = np.zeros((4, N))
sol_4 = np.zeros((4, N))
sol_1[:,0] = c_1
sol_2[:,0] = c_2
sol_3[:,0] = c_3
sol_4[:,0] = c_4

#on boucle sur chaque pas de temps
def resolution(sol):
    for i in range (N-1):
        K1 = h * derivative_1(sol[:,i])
        K2 = h * derivative_1(sol[:,i] + K1/ 2)
        K3 = h * derivative_1(sol[:,i] + K2/ 2)
        K4 = h * derivative_1(sol[:,i]+K3)
        sol[:,i+1] = sol[:,i] + (K1 + (2 * K2) + (2 * K3) + K4)/ 6
        for j in range(2): 
            for k in range(25):
                if abs(sol[j,i+1]) > 2 * (k + 1) * (np.pi):
                    if sol[j,i+1] < 0:
                        sol[j,i+1] += 2 * (k + 1) * (np.pi)
                    else:
                        sol[j,i+1] -= 2 * (k + 1) * (np.pi)
                        break
                else:
                    pass
resolution(sol_1) 
resolution(sol_2)  
resolution(sol_3)  
resolution(sol_4)  

fig = plt.subplots()
plt.plot(t, sol_4[0,:])
plt.title(' oscillation theta c4', fontsize = 25)
plt.xlabel('temps (s)', fontsize = 25)
plt.ylabel('angle (rad)', fontsize = 25)
plt.xlim(35, 45)

fig = plt.subplots()
plt.plot(t, sol_4[1,:])
plt.title(' oscillation phi c4', fontsize = 25)
plt.xlabel('temps (s)', fontsize = 25)
plt.ylabel('angle (rad)', fontsize = 25)
plt.xlim(35, 45)

#matrice des differences entre solutions pour differentes conditions initiales
diff_1 = sol_1 - sol_2
diff_2 = sol_3 - sol_4
for j in range(2):
    for i in range(N-1): 
        for k in range(5):
            if abs(diff_2[j,i+1]) > 2 * (k + 1) * (np.pi):
                if diff_2[j,i+1] < 0:
                    diff_2[j,i+1] += 2 * (k + 1) * (np.pi)
                else:
                     diff_2[j,i+1] -= 2 * (k + 1) * (np.pi)
                     break
            else:
                pass

fig = plt.subplots()
plt.plot(t, diff_1[1,:])
plt.title(' difference angle c1/c2 (phi)', fontsize = 25)
plt.xlabel('temps (s)', fontsize = 25)
plt.ylabel('angle (rad)', fontsize = 25)
plt.xlim(20, 30)

fig = plt.subplots()
plt.plot(t, diff_2[1,:], 'tab:orange')
plt.title(' difference angle c3/c4 (phi)', fontsize = 25)
plt.xlabel('temps (s)', fontsize = 25)
plt.ylabel('angle (rad)', fontsize = 25)
plt.xlim(20, 30)

#Animation
def visualisation(t, vec_1, vec_2):
    fig, ax = plt.subplots()
    x1 = l * np.sin(vec_1[0,:])
    y1 = -l * np.cos(vec_1[0,:])
    x2 = x1 + l * np.sin(vec_1[1,:])
    y2 = y1 - l * np.cos(vec_1[1,:])
    ax.set_xlim(-2 * l, 2 * l)
    ax.set_ylim(-2 * l, 2 * l)
    line1, = plt.plot([0, x1[0], x2[0]], [0, y1[0], y2[0]], linewidth = 2, marker="o", color="blue")
    
    a1 = l * np.sin(vec_2[0,:])
    b1 = -l * np.cos(vec_2[0,:])
    a2 = a1 + l * np.sin(vec_2[1,:])
    b2 = b1 - l * np.cos(vec_2[1,:])
    line2, = plt.plot([0, a1[0], a2[0]], [0, b1[0], b2[0]], linewidth = 2, marker="o", color="red")
    
    def maj(i):
        line1.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
        line2.set_data([0, a1[i], a2[i]], [0, b1[i], b2[i]])
        return line1, line2, 
    
    ani = FuncAnimation(fig, maj, frames=np.arange(0, N, 10), blit=True, interval = 0.5)
    return ani
ani=(visualisation(t, sol_1, sol_2))
HTML(ani.to_html5_video())

#calcul de la jacobienne
def Jacobian(state):
    a, b, c, d = state
    sin = np.sin(a - b)
    cos = np.cos(a - b)
    df_da = ((1)/(1 + sin**2)**2) * ((((-1) * (omega**2)) * (np.sin(b)) * sin * (2 + cos**2)) + ((c**2) * (2 - 3 * cos**2)) - ((d**2) * cos**3) - (2 * omega**2 * np.cos(a) * (1 + sin**2)) + (4 * omega**2 * np.sin(a) * sin * cos))
    df_db = ((1)/(1 + sin**2)**2) * ((((1) * (omega**2)) * (np.sin(b)) * sin * (2 + cos**2)) - ((c**2) * (2 - 3 * cos**2)) + ((d**2) * cos**3) + (omega**2 * np.cos(b) * cos * (1 + sin**2)) - (4 * omega**2 * np.sin(a) * sin * cos))
    df_dc = ((-2) * c * cos * sin)/(1 + sin**2)
    df_dd = ((-2) * d * sin)/(1 + sin**2)
    dg_da = ((1)/(1 + sin**2)**2) * ((((-2) * (omega**2)) * np.sin(a) * sin * (2 + cos**2)) + (d**2 * (1 - 3 * sin**2)) + (2 * c**2 * cos**3) + (2 * (omega**2) * np.cos(a) * cos * (1 + sin**2)) + (4 * (omega**2) * np.sin(b) * sin * cos))
    dg_db = ((1)/(1 + sin**2)**2) * ((((-2) * (omega**2)) * np.sin(a) * sin * (2 + cos**2)) + (d**2 * (1 - 3 * sin**2)) + (2 * c**2 * cos**3) + (2 * (omega**2) * np.cos(a) * cos * (1 + sin**2)) + (4 * (omega**2) * np.sin(b) * sin * cos))
    dg_dc = (4 * c * sin)/(1 + sin**2)
    dg_dd = (2 * d * cos * sin)/(1 + sin**2)
    M = np.array([[0,0,1,0], [0,0,0,1], [df_da,df_db,df_dc,df_dd], [dg_da,dg_db,dg_dc,dg_dd]])
    return M

#Coefficient de Lyapunov
#deviation pour sol_3, cad c3 et sa perturbation
dx = np.zeros((4, N))
dx_0 = [0, (np.pi)/8, 0, 0]
dx[:,0] = dx_0
norme = np.zeros((1, N))
norme[0,:] = np.linalg.norm(dx_0)

def derivative_2(u, v):
    w = np.dot(Jacobian(u), v)
    return w
L = 0
for i in range(N-1):
    dx[:,i+1] = dx[:,i] + h * derivative_2(sol_3[:,i], dx[:,i])
    norme[0,i+1] = np.linalg.norm(dx[:,i+1])
    if i%100 == 0:
        L += np.log(norme[0,i])
lam_max = np.exp(L)

         
            
        





    















    
    
    

    
    
    







    
    


    
    
    

    
    
    








    
    
    

    
    
    






