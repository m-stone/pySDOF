# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:27:14 2021

@author: stone
"""

import numpy as np
import matplotlib.pyplot as plt

def NewmarkBeta(Fext, Pn, m, k, omega, xi, gamma, beta, delta_t):
    # TODO -- non-zero initial conditions, inelastic response (12/1)
    # create output matrices u,v,a
    u,v,a = np.zeros((3,Fext.size+1))
    # initial acceleration
    a[0] = Fext[0] / m - 2*omega*xi*v[0] - omega**2 * u[0]
    
    # Inelastic response
    # determine elastic limits
    u_elastic = Pn / k
    uy_pos = u_elastic
    uy_neg = -1*u_elastic 
    # current plastic deformation
    u_plastic = 0
    # Allocate internal force matrix
    Fint = np.zeros(Fext.size)
    
    # General solver
    for idx,force in enumerate(Fext):
        # Calculate v
        v[idx+1] = v[idx] + (1-gamma)*delta_t*a[idx] + gamma*delta_t*a[idx+1]
        # calculate u
        u[idx+1] = u[idx] + delta_t * v[idx] + 0.5*(delta_t**2)*(
            (1-2*beta)*a[idx]+2*beta*a[idx+1])
        # calculate a
        if idx == Fext.size-1:
            a[idx+1] = a[idx]
        else:
            a[idx+1] = Fext[idx+1] / m - 2*omega*xi*v[idx+1] - omega**2 * u[idx+1]
    return [u,v,a]
    
# Test case
# F_ext = np.zeros(500)

# for idx in range(F_ext.size):
#     if idx < int(F_ext.size/20):
#         F_ext[idx] += 0.5
# [u,v,a] = NewmarkBeta(F_ext,0,10,40,2,0.05,0.5,1/4,0.01)
# fig, ax = plt.subplots(4,1,constrained_layout = True)
# ax[0].plot(F_ext)
# ax[1].plot(u)
# ax[2].plot(v)
# ax[3].plot(a)