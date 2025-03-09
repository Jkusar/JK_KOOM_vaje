# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 23:18:32 2025

@author: Jernej Kušar
"""

import numpy as np

M = 1.3 #W/m^2
f_cl = 1.33
I_cl = 0.18 #m^2K/W
p_z = 13.2 * 10**2 #Pa
t_z = 23 #°C
v_z = 1 #m/s
t_g = 20 #°C
t_r = ((t_g + 273)**4 + (2.5 * 10**8) * (v_z**0.6) * (t_g - t_z))**(1/4) - 273 #°C




def bilanca(M, W, f_cl, t_cl, t_r, a_k, t_z, p_z):
    a = ((t_cl + 273)**4) - ((t_r + 273)**4)
    b = 5733 - 6.99 * (M - W) - p_z
    c = (M - W) - 58.15
    d = 0.0014 * M * (34 - t_z)
    e = 1.73*10**(-5) * M * (5867 - p_z)
    return -(3.97*10**(-8) * f_cl * a + a_k * f_cl * (t_cl - t_z) + (3.05 * 10**(-3)) * b + 0.42 * c + d + e) + M
    
    
def temp_obl(M, W, I_cl, f_cl, t_cl, t_r, a_k, t_z):
    a1 = ((t_cl + 273)**4) - ((t_r + 273)**4)
    a = 3.96 * 10**(-8) * f_cl * a1
    return 35.7 - 0.0257 * (M - W) - I_cl * (a + f_cl * a_k * (t_cl - t_z))


#Začetni približek t_cl
t_cl_0 = 23 #°C
W_0 = 100 #W/m^2

#Iteracijski podatki
max_iter = 100000
min_resid = 0.001

it = 0
    
residuals = {"W":[], "it":[]}

while it < max_iter:
    it += 1
    left = 2.38*(t_cl_0-t_z)**(1/4)
    right = 12.1*np.sqrt(v_z)
    
    if left > right:
        a_k = left
    else:
        a_k = right
        

    W_new = bilanca(M, W_0, f_cl, t_cl_0, t_r, a_k, t_z, p_z)
    
    
    W_resid = np.abs(W_0 - W_new)
    
    residuals["W"].append(W_resid)
    residuals["it"].append(it)
    
    W_0 = W_new

    if W_resid <= min_resid:
        break
    
t_cl_new = temp_obl(M, W_new, I_cl, f_cl, t_cl_0, t_r, a_k, t_z)

import matplotlib.pyplot as plt 
plt.rcParams['figure.dpi'] = 1000

j = -1

plt.plot(residuals["it"][:j], residuals["W"][:j], "b")
plt.xlabel("Iterations")
plt.ylabel("Residuals")
plt.grid()
plt.show


print(f"Efektivna mehanska moč: {np.round(W_new, 4)} W/m^2")
print(f"Ostanek: {np.round(W_resid, 5)}")
print(f"Površinska temperatura obleke: {np.round(t_cl_new, 4)} °C")

