# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 13:20:48 2021

@author: stone
"""

import numpy as np

def CDMnonlin(Fext, Pn, h, omega, xi, m, t_stop = -1):
    # check end time if t_end != -1
    if t_stop != -1:
        # locate stop index
        i_stop = int(t_stop / h)
    else:
        i_stop = Fext.size + 1
    # predefine constants
    k = omega**2 * m
    # create u,v,a matrixes
    u,v,a = np.zeros((3,Fext.size+1))
    # reset all non-inital values to NaN
    u[1:],v[1:],a[1:] = np.nan,np.nan,np.nan
    # assign initital acceleration
    a[0] = Fext[0] / m
    # determine elastic limits
    u_elastic = Pn / k
    uy_pos = u_elastic
    uy_neg = -1*u_elastic
    # current plastic deformation
    u_plastic = 0
    # Allocate internal force matrix
    Fint = np.zeros(Fext.size)
    # reset to nan
    Fint[:] = np.nan
    u_previous = u[0] - (h*v[0]) + (0.5 * (h**2) * a[0])
    for idx, force in enumerate(Fext):
        # check to see if we should stop
        if idx == i_stop:
            break
        if u[idx] > uy_pos:
            # plastic loading in positive direction
            Fint[idx] = Pn
            uy_pos = u[idx]
            u_plastic = uy_pos - u_elastic
            uy_neg = u_plastic - u_elastic
        elif u[idx] < uy_neg:
            # plastic loading in negative direction
            Fint[idx] = -1*Pn
            uy_neg = u[idx]
            u_plastic = uy_neg + u_elastic
            uy_pos = u_plastic + u_elastic
        else:
            # elastic loading
            Fint[idx] = (u[idx] - u_plastic) * k
        # calculate u_previous
        if idx == 0:
            pass
        else:
            u_previous = u[idx-1]
        # calculate u,v,a_i+1:
        u[idx+1] = ( force/m - Fint[idx]/m + (2/h**2)*u[idx] - 
                    (1/h**2 - xi*omega/h)*u_previous ) / (1/h**2 + xi*omega/h)
        v[idx+1] = (u[idx+1] - u_previous) / (2*h)
        a[idx+1] = (u_previous - 2*u[idx] + u[idx+1])/(h**2)
    return [u,v,a,Fint,u_plastic]
    

def CDMnonlinInitValue(Fext,Pn,h,omega,xi,m,
                       u0,v0,a0,Fint,
                       t_start = -1, u_plastic = 0):
    """
    

    Parameters
    ----------
    Fext : NDARRAY
        Input array of external force.
    Pn : float
        Yield force.
    h : float
        step size of arrays.
    omega : float
        structural frequency.
    xi : float
        structural damping factor.
    m : float
        structural mass.
    u0 : ndarray
        displacement array.
    v0 : ndarray
        velocity array.
    a0 : ndarray
        acceleration array.
    Fint : ndarray
        internal force array.
    t_start : float, optional
        re-start time. The default is -1.
    u_plastic : float, optional
        initial plastic deformation. The default is 0.

    Returns
    -------
    list
        Returns list of 3 ndarrays and 1 float:
            [u, v, a, F_int, u_plastic]

    """
    print("CDM Nonlin Initial Values (Restart Analysis)")
    if t_start != -1:
        i_start = int(t_start/h)
    else:
        i_start = 0
    # predefine constants
    k = (omega**2) * m
    # create u,v,a matrixes
    u,v,a = u0, v0, a0
    # determine elastic limits
    u_elastic = Pn / k
    # if internal force is past the elastic limit:
    if np.abs(Fint[i_start]) >= Pn:
        if np.abs(Fint[i_start]) >= Pn:
            uy_pos = u[i_start+1]
            u_plastic = uy_pos - u_elastic
            uy_neg = u_plastic - u_elastic
        else:
            uy_neg = u[i_start+1]
            u_plastic = uy_neg + u_elastic
            uy_pos = u_plastic + u_elastic
    else:
        # we are in elastic range, so go ahead and reset the uy limits
        uy_pos = u_plastic + u_elastic
        uy_neg = u_plastic - u_elastic
    # loop through external force
    for idx,force in enumerate(Fext):
        # if we are before start time, do nothing
        if idx < i_start:
            pass
        # else actually calculate changes to u,v,a
        else:
            if u[idx] > uy_pos:
                # plastic loading in positive direction
                Fint[idx] = Pn
                uy_pos = u[idx]
                u_plastic = uy_pos - u_elastic
                uy_neg = u_plastic - u_elastic
            elif u[idx] < uy_neg:
                # plastic loading in negative direction
                Fint[idx] = -1*Pn
                uy_neg = u[idx]
                u_plastic = uy_neg + u_elastic
                uy_pos = u_plastic + u_elastic
            else:
                # elastic loading
                Fint[idx] = (u[idx] - u_plastic) * k
            # calculate u_previous
            if idx == 0:
                pass
            else:
                u_previous = u[idx-1]
            # calculate u,v,a_i+1:
            u[idx+1] = ( force/m - Fint[idx]/m + (2/h**2)*u[idx] - 
                        (1/h**2 - xi*omega/h)*u_previous ) / (1/h**2 + xi*omega/h)
            v[idx+1] = (u[idx+1] - u_previous) / (2*h)
            a[idx+1] = (u_previous - 2*u[idx] + u[idx+1])/(h**2)
    return [u,v,a,Fint,u_plastic]

"""
Tedesco ex15.1

# structural properties
# given properties
k = 12.35
c = 0.27
m = 0.20
Rm = 15.0

# derived properties
# w = sqrt(k/m)
omega = np.sqrt(k/m)
# xi = c/c_crit = c / (2*m*w)
xi = c/(2*m*omega)

# solution properties
t_f = 3.0
h = 0.001
n = int(t_f/h)
t = np.arange(0.0,t_f,h)
t_plot = np.append(t,t_f+h)

# external force
f_t = np.array([0.0,20.0,0.0,-10.0,0.0,0.0])
f_x = np.array([0.0,0.45,1.1,1.2,1.4,t_f])
F_ext = np.interp(t,f_x,f_t)

[u,v,a] = np.zeros((3,n+1))
f_int = np.zeros(n)
"""

   
