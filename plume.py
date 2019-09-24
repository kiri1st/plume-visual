##############################################################
## FUNCTIONS FOR CALCULATING PLUME VARIABLES 
##############################################################
import math 
import numpy as np
##
## 
## Fxn calculates u*, the friction velocity
def u_star(M_10, z_R=10, z_0=1.0):
    ## Equation 18.13
    ## Friction velocity in terms of surface wind speed and roughness length 
    # See Table 18-1 for z_0 values:
    # z_0 = 0.1 # [m]  this is for 'roughly open terrain'
    # z_0 = 1.0 # [m] this is for 'closed terrain', suburbia

    ## In statistically neutral conditions:
    k = 0.4 # von Karman constant
    # z_R = 10 # [m] standard height for measuring surface winds
    # M_10 = 10 # [m/s] wind speed at 10 m
    u_star = (k * M_10)/(np.log(z_R/z_0))  # friction velocity
    
    return u_star
##
##
## Fxn returns sigma_v and sigma_w, the friction velocity
def turbulence(z_CL, h, M_10, z_R=10, z_0=1.0):
    ## Variables:
        # z = height of detector ## <- nope - Stull says height of centerline
        # h = depth of ABL
    ## Variables that go to u_star():
        # M_10 = windspeed near surface
        # z_R = height of near surface windspeed
        # z_0 = friction coefficient of surface
    
    ## Equation 18.25 b & c
    # for neutral stablity:
    u = u_star(M_10, z_R=z_R, z_0=z_0)
    sigma_v = 1.6 * u * (1 - (0.5 * (z_CL / h)))
    sigma_w = 1.25 * u * (1 - (0.5 * (z_CL / h)))
    
    return sigma_v, sigma_w
##
##
## Fxn calculates the centerline height of a plume
# NOTE: BE CAREFUL IF YOU USE THIS - VERIFY THAT IT WORKS
#       and be sure to check other functions that depend on z_CL
def centerline_height(x, z_s=10, R_0=0.5, W_0=10, M=5, 
                      theta_a=270, theta_p=300):
    # variables:
        # R_0 = stack diameter [m]
        # W_0 = exit velocity of effluent [m/s]
        # M = ambient wind speed at stack top [m/s]
        # z_s = stack height [m]
        # theta_a = ambient potential temp
        # theta_p = potential temp of plume
    
    # constants
    a = 8.3; b = 4.2
    
    # momentum length scale
    l_m = W_0 * R_0 / M 
    
    # temperature excess of effluent
    delta_theta = (theta_p - theta_a)/theta_a
    
    # buoyancy length scale (l_b)
    l_b = ((W_0 * R_0**2 * 9.8) / M**3)*delta_theta
    
    # solve for the plume centerline height
    z_CL = z_s + ((a * l_m**2 * x) + (b * l_b * x**2))**(1/3)
#     z_CL = z_s # for now, z_CL is just the height of the stack
    
    return z_CL
##
##
## Fxn calculates dispersion of a plume in the y direction
def y_disp(U, x, y, sigma_v=0.96):
    ## Check conditions to avoid division by zero
    if x == 0: 
        x = 0.1 # make x close to zero to avoid an error
#         sigma_y = 1

        ## Calculate sigma_y
        t_L = 60 # Lagrangian time scale [sec]
        # Stull eq 19.14 (close approximation)
    #     sigma_y = sigma_v * x / U  # [m]

        # Stull eq 19.15 (far approximation)
    #     sigma_y = sigma_v * (2 * t_L * (x / U))**(1/2)  # [m]

        # Stull eq 19.13a
        sigma_y = math.sqrt(2 * sigma_v**2 * t_L **2 *((x / (U*t_L)) - \
                        1 + np.exp(-(x / (U*t_L)))))
    else:
        ## Calculate sigma_y
        t_L = 60 # Lagrangian time scale [sec]
        # Stull eq 19.14 (close approximation)
    #     sigma_y = sigma_v * x / U  # [m]

        # Stull eq 19.15 (far approximation)
    #     sigma_y = sigma_v * (2 * t_L * (x / U))**(1/2)  # [m]

        # Stull eq 19.13a
        sigma_y = math.sqrt(2 * sigma_v**2 * t_L **2 *((x / (U*t_L)) - \
                        1 + np.exp(-(x / (U*t_L)))))

    ## Use sigma_y to calculate disp_y
    disp_y = np.exp(-0.5*((y/sigma_y)**2))
    
    return disp_y, sigma_y
##
##
## Fxn calculates dispersion of a plume in the z direction
def z_disp(U, x, z, z_CL = 10, sigma_w=1.04):
    #     z_CL = centerline_height(x)
    
    ## Check conditions to avoid division by zero
    if x == 0:
        sigma_z = 1
       
        ## Use sigma_z to calculate disp_z
        if z == 0:
            disp_z = np.exp(-0.5*(((z_CL)/sigma_z)**2))
        else:
            disp_z = np.exp(-0.5*(((z - z_CL)/sigma_z)**2)) + \
                    np.exp(-0.5*(((z + z_CL)/sigma_z)**2))
    else:
        ## Calculate sigma_z
        t_L = 60 # Lagrangian time scale [sec]
        # Stull eq 19.14
    #     sigma_z = sigma_w * x / U  # [m]

        # Stull eq 19.15
    #     sigma_z = sigma_w * (2 * t_L * (x / U))**(1/2)  # [m]

        # Stull eq 19.13a
        sigma_z = math.sqrt(2 * sigma_w**2 * t_L **2 *((x / (U*t_L)) - \
                         1 + np.exp(-x/(U*t_L))))

        ## Use sigma_z to calculate disp_z
        if z == 0:
            disp_z = np.exp(-0.5*(((z_CL)/sigma_z)**2))
        else:
            disp_z = np.exp(-0.5*(((z - z_CL)/sigma_z)**2)) + \
                    np.exp(-0.5*(((z + z_CL)/sigma_z)**2))

    return disp_z, sigma_z
##
##
## Fxn calculates the concentration at given coordinates
def C_EtO(Q, U, x, y, z, z_CL=10, sigma_v=0.8, sigma_w=1.04):
    ## note: U = M, wind speed at the plume centerline height

    # dispersion in y and z directions
    y_dispersion, sigma_y = y_disp(U, x, y, sigma_v=sigma_v)
    z_dispersion, sigma_z = z_disp(U, x, z, z_CL, sigma_w=sigma_w)

    # advection and dilution with distance x
    if z == 0:
        advection = math.pi*sigma_y*sigma_z*U
    else: 
        advection = 2*math.pi*sigma_y*sigma_z*U

    # solve full equation
    C_EtO = (Q/advection) * y_dispersion * z_dispersion
    
    return C_EtO
##   
##
#############################################################################
## PLOTTING FXNS
#############################################################################
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
##
##
## Fxn preps data set for plotting plume
def plume_2Dmesh(h, z_CL, M_10, Q, U, z=0, z_0=1.0, xmax = 2000, ymax = 500):
    sigma_v, sigma_w = turbulence(z_CL, h, M_10, z_0=z_0)
    
    ## Make 2D arrays representing the surface area using meshgrid
    n = 100
    x_array = np.linspace(0,xmax,num=n)
    y_array = np.linspace(-ymax,ymax,num=n)
    X, Y = np.meshgrid(x_array, y_array)
    C = np.zeros_like(X) 
    for i in range(n):
        for j in range(n):
            conc = C_EtO(Q, U, X[i,j], Y[i,j], z, z_CL = z_CL, sigma_v=sigma_v,                                          sigma_w=sigma_w)
            C[i,j] = conc
    return X, Y, C
##
##
## Fxn plots a 2D concentration plume
def plot_plume(X, Y, C, maxC, fig, ax):
    ## Match a color map to the concentration range
    levels = MaxNLocator(nbins=40).tick_values(0, maxC)
    cmap = plt.get_cmap('magma_r')
#     cmap = plt.get_cmap('YlOrBr') # alternate color scheme
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    
    ## Make the axis object and figure
    cf = ax.contourf(X, Y, C, levels=levels, cmap=cmap)
    fig.colorbar(cf, ax=ax)

    return fig, ax
##