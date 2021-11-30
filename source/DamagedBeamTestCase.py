# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 13:59:25 2021
Damaged beam test case. Use Feldman-Seiss single reinforced, with damage
removing top 2inches of cover.
b = 6in
h = 12 in
d = 10.5in
As = 1.20 in^2
f`c = 5300 psi
@author: stone
"""

import numpy as np
import matplotlib.pyplot as plt
import CDMnonlin
import pandas as pd # pandas arrays to export to DPLOT. 
from pathlib import Path


# undamaged Slab Values
Py_u = 21.0
xy_u = 1.0
K_u = Py_u / xy_u

m_u = 1.20012200011
xi = 0.05

# undamaged Resistance values - from DSAS
P_y_u = 21.0455
x_y_u = 0.243651
K_u = P_y_u / x_y_u

m_u = 0.750 / 386.1 #kips-s^2 / in
xi = 0.05

omega_u = np.sqrt(K_u / m_u)
c_crit_u = 2*m_u *omega_u
c_u = xi * c_crit_u

# damaged resistance values
P_y_d = 15.3577
x_y_d = 0.296935
K_d = P_y_d / x_y_d

m_d = 0.625 / 386.1 #kip-s^2/in
xi = 0.05

omega_d = np.sqrt(K_d / m_d)
c_crit_d = 2*m_d *omega_d
c_d = xi * c_crit_d

# solution properties
t_f = 0.500
h = 0.0001
n = int(t_f/h)
t = np.arange(0.0,t_f,h)
t_plot = np.append(t,t_f+h)

# external force array
# "blast"
f_y = np.array([0.0,25.00,0.000,-5.00,0.000,0.0]) #Load
f_x = np.array([0.0,0.010,0.030,0.035,0.040,t_f]) #Time
# "fragment" -- Note that this is added to f
f2_y = np.array([0.000,030.0,0.000]) # Load
f2_x = np.array([0.013,0.015,0.017]) # time
F_ext = np.interp(t,f_x,f_y) + np.interp(t,f2_x,f2_y)

# Plot load vs time
fig0, ax0 = plt.subplots()
ax0.plot(t[:500],F_ext[:500])
ax0.set_title('Load vs Time')
ax0.set_ylabel('Load')
ax0.set_xlabel('Time')
ax0.grid()

# output variables
# undamaged 
[u_u,v_u,a_u] = np.zeros((3,n+1))
f_int_u = np.zeros(n)

# damaged
[u_d,v_d,a_d] = np.zeros((3,n+1))
u_d[:] = np.nan
f_int_d = np.zeros(n)

# run CDM case (no damage)
print("Running CDM Nonlin with no damage evolution:")
[u_u,v_u,a_u,f_int_u,u_plastic] = CDMnonlin.CDMnonlin(F_ext, P_y_u, h, omega_u, xi, m_u)
print("Final plastic deformation: {}".format(u_plastic))
#[u_d,v_d,a_d,f_int_d] = CDMnonlin.CDMnonlin(F_ext, P_y_d, h, omega_d, xi, m_d)
#fig2, axs = plt.subplots(2,1)
#axs[1].plot(t_plot,u)
fig, ax = plt.subplots(5,1,constrained_layout = True)
fig.set_size_inches(8.5,11)
ax[0].plot(t,F_ext)
ax[1].plot(t,f_int_u)
ax[2].plot(t_plot,u_u)
ax[3].plot(t_plot,v_u)
ax[4].plot(t_plot,a_u)
ax[0].set_title('Load vs Time')
ax[1].set_title('Internal Force vs Time (undamaged)')
ax[2].set_title('Displacement vs Time (undamaged)')
ax[3].set_title('Velocity vs Time (undamaged)')
ax[4].set_title('Acceleration vs Time (undamaged)')
for axis in ax:
    axis.grid()
fig.suptitle('Undamaged CDM Example')


#so want to run CDMNonlin up to t_d, then extract x,v,a(,F_int?). 
print("Running CDM Nonlin with damage evolution:")
# now assume damage occurs at t = 0.017s
t_damage = 0.0150
idx_damage = int(t_damage / h)
print("Damage / restart occurs at t = {}".format(t_damage))
print("Running until restart time. . .")
[u_d,v_d,a_d,f_int_d,u_plastic] = CDMnonlin.CDMnonlin(F_ext, P_y_u, h, omega_u, xi, m_u,t_damage)
print("Finished initial analysis.")
print("Current plastic deformation: {}".format(u_plastic))
print("u: {:.3f}".format(u_d[idx_damage]))
print("v: {:.3f}".format(v_d[idx_damage]))
print("a: {:.3f}".format(a_d[idx_damage]))

print("Restarting with damage. . .")
# then restart analysis with existing initial values

[u_d,v_d,a_d,f_int_d,u_plastic] = CDMnonlin.CDMnonlinInitValue(F_ext, P_y_d, h,
                                                               omega_d, xi, m_d,
                                                               u_d,v_d,a_d,f_int_d,
                                                               t_start = t_damage,
                                                               u_plastic = u_plastic)

print("Final plastic deformation: {}".format(u_plastic))
fig2, ax2 = plt.subplots(5,1,constrained_layout=True)
fig2.set_size_inches(8.5,11)
ax2[0].plot(t,F_ext)
ax2[1].plot(t,f_int_d)
ax2[2].plot(t_plot,u_d)
ax2[3].plot(t_plot,v_d)
ax2[4].plot(t_plot,a_d)
ax2[0].set_title('Load vs Time')
ax2[1].set_title('Internal Force vs Time (damaged)')
ax2[2].set_title('Displacement vs Time (damaged)')
ax2[3].set_title('Velocity vs Time (damaged)')
ax2[4].set_title('Acceleration vs Time (damaged)')
for axis in ax2:
    axis.grid()
fig2.suptitle('Damaged CDM Restart')

# Create Dataframes for DPLOT
# first, resize F_ext and F_int so its on same timebase
# append last value
F_ext = np.append(F_ext,F_ext[-1])
f_int_u = np.append(f_int_u,f_int_u[-1])
f_int_d = np.append(f_int_d, f_int_d[-1])
# create undamaged dataframe
df_undamaged = pd.DataFrame({'Time':t_plot,
                             'F{\dext}':F_ext,
                             'F{\dint}':f_int_u,
                             'u{\dundamaged}':u_u,
                             'v{\dundamaged}':v_u,
                             'a{\dundamaged}':a_u })
# reset index to time column
df_undamaged.set_index('Time',inplace=True)
# create damaged dataframe
df_damaged = pd.DataFrame({'Time':t_plot,
                             'F{\dext}':F_ext,
                             'F{\dint}':f_int_d,
                             'u{\ddamaged}':u_d,
                             'v{\ddamaged}':v_d,
                             'a{\ddamaged}':a_d })
# reset index to time column
df_damaged.set_index('Time',inplace=True)

# output to CSV if necessary, or just copy and paste from Spyder
fname_u=Path('undamaged.csv')
fname_d=Path('damaged.csv')
df_undamaged.to_csv(fname_u)
df_damaged.to_csv(fname_d)

# u_curve = dp.curve(t_plot,u_d,title='u_damaged')
# dplt = dp.plot(u_curve,xlabel='Time(s)',ylabel='Displacement(in)')
# dplt.show()

