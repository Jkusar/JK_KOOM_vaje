# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 23:18:32 2025

@author: Jernej Kušar
"""

import numpy as np

M = 1.3*58.2 #W/m^2
f_cl = 1.33
I_cl = 0.17 #m^2K/W
p_z = 13.2 * 10**2 #Pa
t_z = 23 #°C
v_z = 0.1 #m/s
t_g = 23.5 #°C
t_r = ((t_g + 273)**4 + (2.5 * 10**8) * (v_z**0.6) * (t_g - t_z))**(1/4) - 273 #°C


def bilanca(M, H, f_cl, t_cl, t_r, a_k, t_z, p_z):
    a = ((t_cl + 273)**4) - ((t_r + 273)**4)
    b = 5733 - 6.99 * H - p_z
    c = H - 58.15
    d = 0.0014 * M * (34 - t_z)
    e = 1.73*10**(-5) * M * (5867 - p_z)
    return 3.97*10**(-8) * f_cl * a + a_k * f_cl * (t_cl - t_z) + (3.05 * 10**(-3)) * b + 0.42 * c + d + e
    
    
def temp_obl(M, H, I_cl, f_cl, t_cl, t_r, a_k, t_z):
    a1 = ((t_cl + 273)**4) - ((t_r + 273)**4)
    a = (3.96 * 10**(-8) * f_cl * a1)
    return 35.7 - 0.0275 * H - I_cl * (-a + (f_cl * a_k * (t_cl - t_z)))


def PMV(M, H, f_cl, t_cl, t_r, t_z, a_k, p_z):
    a = (0.303 * np.exp(-0.036 * M) + 0.028)
    b1 = ((t_cl + 273)**4) - ((t_r + 273)**4)
    b = 3.96 * 10**(-8) * f_cl * b1 
    c = a_k * f_cl * (t_cl - t_z)
    d = 3.05 * 10**(-3) * (5733 - 6.99 * H - p_z)
    e = 0.42 * (H - 58.15)
    f = 0.0014 * M * (34 - t_z)
    g = 1.7 * 10**(-5) * M * (5867 - p_z)
    return a*(H-b-c-d-e-f-g)



def PPD(PMV):
    return 100 - 95*np.exp(-0.03353 * (PMV**4) - 0.2179 * PMV**2)



#Začetni približek t_cl
t_cl_0 = 23 #°C
H_0 = M #W/m^2

#Iteracijski podatki
max_iter = 100000
min_resid = 0.001

it = 0
    
residuals = {"H":[], "t":[], "it":[]}

while it < max_iter:
    
    it += 1
    left = 2.38*(t_cl_0-t_z)**(1/4)
    right = 12.1*np.sqrt(v_z)
    
    if left > right:
        a_k = left
    else:
        a_k = right
            
    
        
    H_new = bilanca(M, H_0, f_cl, t_cl_0, t_r, a_k, t_z, p_z)
    t_cl_new = temp_obl(M, H_new, I_cl, f_cl, t_cl_0, t_r, a_k, t_z)

    
    
    H_resid = np.abs(H_0 - H_new)
    t_resid = np.abs(t_cl_0-t_cl_new)
    
    residuals["H"].append(H_resid)
    residuals["t"].append(t_resid)
    residuals["it"].append(it)
    
    H_0 = H_new
    t_cl_0 = t_cl_new

    if t_resid <= min_resid and H_resid <= min_resid:
        break



import matplotlib.pyplot as plt 
plt.rcParams['figure.dpi'] = 1000

j = -1

plt.plot(residuals["it"][:j], residuals["t"][:j], "b", label = "T")
plt.plot(residuals["it"][:j], residuals["H"][:j], "orange", label = "W")
plt.xlabel("Iterations")
plt.ylabel("Residuals")
plt.legend()
plt.grid()
plt.show


print(f"Efektivna mehanska moč: {np.round(H_new, 4)} W/m^2")
print(f"Ostanek H: {np.round(H_resid, 5)}")

print(f"Površinska temperatura obleke: {np.round(t_cl_new, 4)} °C")
print(f"Ostanek t: {np.round(t_resid, 5)}")



PMV_res = PMV(M, H_new, f_cl, t_cl_new, t_r, t_z, a_k, p_z)

PPD_res = PPD(PMV_res)


print(f"PMV: {np.round(PMV_res, 4)}")
print(f"PPD: {np.round(PPD_res, 4)}")



