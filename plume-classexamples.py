###############################################
#### PLUME PLOTTING EXAMPLES FROM CLASS ####
###############################################

## import everything from plume.py
from plume import *
# imports math, numpy as np, matplotlib.pyplot as plt
# imports functions that calculate plume variables 
# imports functions for plotting: plume_2Dmesh and plot_plume 

###############################################
## Effect of Q on plume concentration
###############################################
# define variables:
h = 1500 # approx depth of ABL 
z_CL = 10 # height of centerline - assume constant
U = 4 # windspeed at centerline
M_10 = U # since centerline is also 10m, these should be the same
z = 4 # approx height of monitoring stations

# pick two emission velocities
Q1 = 100000 # ug/L
Q2 = 200000 # ug/L

# get the 2D arrays
X1, Y1, C1 = plume_2Dmesh(h, z_CL, M_10, Q1, U, z=z, xmax = 500, ymax = 100)
X2, Y2, C2 = plume_2Dmesh(h, z_CL, M_10, Q2, U, z=z, xmax = 500, ymax = 100)

# get the max concentration from both sets so plots have the same legend
if C2.max() > C1.max():
    maxC = C2.max()
else: 
    maxC = C1.max()

# make the plume figures 
fig, (ax0, ax1) = plt.subplots(nrows=2,figsize=(10,10))

# plot plume 1
fig, ax0 = plot_plume(X1, Y1, C1, maxC, fig, ax0)
ax0.set_title('Plume with Q = 100,000 ug/s')
ax0.set_ylabel('Distance y from centerline x (m)')

# plot plume 2
fig, ax1 = plot_plume(X2, Y2, C2, maxC, fig, ax1)
ax1.set_title('Plume with Q = 200,000 ug/s')
ax1.set_xlabel('X direction - distance from source (m)')
ax1.set_ylabel('Distance y from centerline x (m)')

fig.tight_layout()
plt.show()

################################################
## Effect of wind speed on plume concentration
################################################
# define variables:
h = 1500 # approx depth of ABL 
z_CL = 10 # height of centerline - assume constant
z = 4 # approx height of monitoring stations
Q = 200000 # ug/L 

# pick two wind speeds
U1 = 4
U2 = 11

# get the 2D arrays
X1, Y1, C1 = plume_2Dmesh(h, z_CL, U1, Q, U1, z=z, xmax = 500, ymax = 100)
X2, Y2, C2 = plume_2Dmesh(h, z_CL, U2, Q, U2, z=z, xmax = 500, ymax = 100)

# get the max concentration from both sets so plots have the same legend
if C2.max() > C1.max():
    maxC = C2.max()
else: 
    maxC = C1.max()
    
# make the plume figures  
fig, (ax0, ax1) = plt.subplots(nrows=2,figsize=(10,10))

# plot plume 1
fig, ax0 = plot_plume(X1, Y1, C1, maxC, fig, ax0)
ax0.set_title('Plume with U = 4 m/s')
ax0.set_ylabel('Distance y from centerline x (m)')

# plot plume 2
fig, ax1 = plot_plume(X2, Y2, C2, maxC, fig, ax1)
ax1.set_title('Plume with U = 11 m/s')
ax1.set_xlabel('X direction - distance from source (m)')
ax1.set_ylabel('Distance y from centerline x (m)')

fig.tight_layout()
plt.show()

#############################################################
## Effect of wind speed - vertical concentration profile
#############################################################
# define variables:
h = 1500 # approx depth of ABL 
z_CL = 20 # height of centerline - assume constant
U = 9 # windspeed at centerline
M_10 = U # since centerline is also 10m, these should be the same
Q = 20000 # ug/L # avg emission times number of stacks  
x = 100
y = 0

# initialize figure:
fig = plt.figure(figsize=(10,5))
plt.xlabel('concentration (ug/L)', fontsize=12)
plt.ylabel('height z (m)', fontsize=12)
plt.title('Change in vertical concentration profile with wind speed')

# pick a range of wind speeds (U)
colors = ['r','y','g','c','b','k','m']; c = 0
for U in np.linspace(1,15,7):
    # plot the vertical concentration profile for each wind speed
    for z in np.linspace(1,40):
        sigma_v, sigma_w = turbulence(z_CL, h, U)
        Cz = C_EtO(Q, U, x, y, z, z_CL=z_CL, sigma_v=sigma_v, sigma_w=sigma_w)
        plt.scatter(Cz, z, c = colors[c])
    print("Windspeed %.1f m/s in color %s" % (U, colors[c]))
    c += 1

