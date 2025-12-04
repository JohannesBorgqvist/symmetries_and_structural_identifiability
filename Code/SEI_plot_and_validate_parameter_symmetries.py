# =========================================================
# SEI validation
# Date: 2024-02-18,
# Written by: Johannes
# =========================================================
# LOAD LIBRARIES
import numpy as np  # For numerical calculations
import matplotlib.pyplot as plt  # For plotting
from scipy.integrate import odeint  # For solving ODEs numerically
# Make the plotting fancy with LaTeX
plt.rcParams['text.usetex'] = True
# =========================================================
# Functions
# =========================================================
#----------------------------------------------------------
# Function 1: Solving our beloved SEI ODE system for the
# putputs
#----------------------------------------------------------
def SEI_outputs(y_vec, t,p):
    # Unpack output states
    x = y_vec[0] # Corresponds to y_I
    y = y_vec[1] # Corresponds to y_E
    z = y_vec[2] # Corresponds to \dot{y_I}
    # Unpack parameters
    c = p[0]
    beta = p[1]
    delta = p[2]
    eps = p[3]
    muS = p[4]
    muE = p[5]
    muI = p[6]
    kE = p[7]
    kI = p[8]
    # Define our ODEs
    dydt = []
    # ODE for \dot{x}
    dydt.append(z)
    # ODE for \dot{y}
    dydt.append(-delta*y/eps - kE*muI*x/kI - kE*z/kI - muE*y + kE*muI*x/(eps*kI) + kE*z/(eps*kI))
    # ODE for \dot{z}
    dydt.append(beta*c*eps*x + beta*delta*x*y/kE - beta*muI*x**2/kI - beta*x*z/kI - delta**2*kI*y/(eps*kE) - delta*muI*x - delta*z - delta*kI*muE*y/kE + delta*kI*muS*y/kE - delta*kI*y*z/(kE*x) + delta*muI*x/eps + delta*z/eps - muI*muS*x - muS*z + z**2/x)
    # Return the ODE system
    return dydt
#----------------------------------------------------------
# Function 2: Parameter symmetry
#----------------------------------------------------------
def param_symmetry_SEI(y, var_eps, alpha):
    # Unpack parameters
    c = y[0]
    beta = y[1]
    delta = y[2]
    eps = y[3]
    muS = y[4]
    muE = y[5]
    muI = y[6]
    kE = y[7]
    kI = y[8]
    # Unpack our two alphas as well
    alpha_1 = alpha[0]
    alpha_2 = alpha[1]
    # Allocate our ODE system for the parameters
    dtheta_dvareps = []
    # Add ODEs to our parameter ODE system
    dtheta_dvareps.append((alpha_2-alpha_1)*c*eps-alpha_2*c) # ODE for c
    dtheta_dvareps.append(alpha_1*beta) # ODE for beta    
    dtheta_dvareps.append((alpha_2-alpha_1)*delta) # ODE for delta    
    dtheta_dvareps.append((alpha_1-alpha_2)*eps*(eps-1)) # ODE for eps    
    dtheta_dvareps.append(0) # ODE for muS
    dtheta_dvareps.append(-(alpha_2-alpha_1)*delta) # ODE form muE    
    dtheta_dvareps.append(0) # ODE for muI
    dtheta_dvareps.append(alpha_2*kE) # ODE for kE    
    dtheta_dvareps.append(alpha_1*kI) # ODE for kI        
    # Return all ODE system for the parameters
    return dtheta_dvareps
# =========================================================
# Solve ODE system
# =========================================================
# State the number of data point
num_data_points = 200
# Simulate time forwards
t = np.linspace(0, 100, num_data_points)
# Define parameters
p = [200, 1.0, 0.1, 0.05, 0.005, 0.005, 0.01, 0.5, 0.5]
# Define alphas for the parameter symmetries
alpha = [2, 1]
# Define ICS of the original SEI model
S0 = 10
E0 = 1
I0 = 1
# Define the ICS
y0 = [E0*p[-2], I0*p[-1], p[-1]*(p[2]*p[1]*S0*I0+p[3]*E0-p[6]*I0)]
# Define the ODEs
SEI_states = odeint(SEI_outputs, y0, t, args=(p,))
# Extract states
yI, yE, dot_yI = SEI_states.T
# ================================================================================================================================================
# Validation of parameter symmetries
# ================================================================================================================================================
# Define a var_eps vector
var_eps_vec = np.linspace(0,1,100)
# Allocate memory for our SSs
SS_yE = []
SS_yI = []
for var_eps in var_eps_vec:
    # Define a var_eps vector
    var_eps_vec_temp = np.linspace(0,var_eps,1000000)
    # Solve ODE system for the parameters
    p_hat_sys = odeint(param_symmetry_SEI, p, var_eps_vec_temp, args=(alpha,))
    # Update the transformed parameter vector
    p_hat = p_hat_sys[-1,0:len(p)]
    # Define the ODEs
    SEI_states = odeint(SEI_outputs, y0, t, args=(p_hat,))
    # Extract states
    yI_hat, yE_hat, dot_yI_hat = SEI_states.T
    # Calculate residuals
    SS_yE.append(np.sqrt(np.inner(yE_hat-yE,yE_hat-yE)/len(yI)))
    SS_yI.append(np.sqrt(np.inner(yI_hat-yI,yI_hat-yI)/len(yI)))
# Re-cast lists for the outputs yI and yE as arrays
SS_yE = np.array(SS_yE)
SS_yI = np.array(SS_yI)
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Plot SEI solutions
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#Define the first figure
f_SEI, ax_SEI = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Subplot A
#---------------------------------------------------------------------------------
# Plot the simulated data and the underlying model
ax_SEI[0].plot(t,yE,color=(0/256,68/256,27/256),label="$y_{E}(t)$",linewidth=3.0)
ax_SEI[0].plot(t,yI,color=(35/256,139/256,69/256),label="$y_{I}(t)$",linewidth=3.0)
# Set the y-limit
#ax_SEI[0].set_ylim([0, 2.5])
# Set a grid and define a legend
ax_SEI[0].grid()
ax_SEI[0].legend(loc='best',prop={"size":20})
# Set the x-labels and y-labels
ax_SEI[0].set_xlabel(xlabel="Time, $t$",fontsize=20)
ax_SEI[0].set_ylabel(ylabel="Population density",fontsize=20)
ax_SEI[0].set_title(label="Dynamics of outputs of SEI-model",fontsize=20)
ax_SEI[0].xaxis.set_tick_params(labelsize=15)
ax_SEI[0].yaxis.set_tick_params(labelsize=15)
#---------------------------------------------------------------------------------
# Subplot B
#---------------------------------------------------------------------------------
# Plot the simulated data and the underlying model
ax_SEI[1].plot(var_eps_vec,SS_yE,color=(2/256,56/256,88/256),label="$\\sqrt{\\frac{\\sum_{i=1}^{n}\\left(y_{E}(t_{i},\\vec{\\theta})-y_{E}(t_{i},\\hat{\\vec{\\theta}}(\\varepsilon))\\right)^{2}}{n}}$",linewidth=3.0)
ax_SEI[1].plot(var_eps_vec,SS_yI,color=(35/256,139/256,69/256),label="$\\sqrt{\\frac{\\sum_{i=1}^{n}\\left(y_{I}(t_{i},\\vec{\\theta})-y_{I}(t_{i},\\hat{\\vec{\\theta}}(\\varepsilon))\\right)^{2}}{n}}$",linewidth=3.0)
# # Set a grid and define a legend
ax_SEI[1].grid()
ax_SEI[1].legend(loc='best',prop={"size":20})
# Set the x-labels and y-labels
ax_SEI[1].set_xlabel(xlabel="$\\varepsilon$",fontsize=20)
ax_SEI[1].set_ylabel(ylabel="Sum of squares, $\\mathrm{SS}(\\varepsilon)$",fontsize=20)
ax_SEI[1].set_title(label="The fit is invariant under transformation by $\\Gamma_{\\varepsilon}^{\\vec{\\theta}}$",fontsize=20)
ax_SEI[1].xaxis.set_tick_params(labelsize=15)
ax_SEI[1].yaxis.set_tick_params(labelsize=15)
#plt.show()
# Title and saving the figure
f_SEI.savefig('../Figures/SEI_numerical_validation.png')

# ================================================================================================================================================
# Plot of action of parameter symmetry (alpha_1,alpha_2)=(2,1) for different start guesses in parameter space
# ================================================================================================================================================
# Have some fancy colours ready in store
clr_vec = [(77/256, 0/256, 75/256), (129/256, 15/256, 124/256), (136/256, 65/256, 157/256), (140/256, 107/256, 177/256), (140/256, 150/256, 198/256), (158/256, 188/256, 218/256), (191/256, 211/256, 230/256), (224/256, 236/256, 244/256), (247/256, 252/256, 253/256)]
# Create a vector p_vec with different starting points in the parameter space
p_vec = []
for index in range(len(clr_vec)):
    p_vec.append(np.array(np.multiply(np.random.rand(1,len(p)),p)[0]))
# Do a var_eps vec for the transformation parameters
var_eps_vec = np.linspace(0,5,100)
# Allocate the number of parameters we want to explore
p_hat = []
# Loop over all alphas and calculate new parameters
for index, p_temp in enumerate(p_vec):
    # Solve ODE system for the parameters
    p_hat_sys = odeint(param_symmetry_SEI, p_temp, var_eps_vec, args=(alpha,))
    # Update the transformed parameter vector
    p_hat.append(p_hat_sys)

#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# 2D parameter plots
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#Define the first figure
f_SEI_2D, ax_SEI_2D = plt.subplots(1, 3, constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------
# Subplot A: delta vs uE
#---------------------------------------------------------------------------------    
for index in range(len(clr_vec)):
    # Plot the parameters against eachother yeah?
    ax_SEI_2D[0].plot(p_hat[index][:,2],p_hat[index][:,5],color=clr_vec[index],linewidth=3.0, linestyle = "dashdot")        
# Set the y-limit
#ax_SEI[0].set_ylim([0, 2.5])
# Set a grid and define a legend
ax_SEI_2D[0].grid()
#ax_SEI_2D[0].legend(loc='best',prop={"size":20})
# Set the x-labels and y-labels
ax_SEI_2D[0].set_xlabel(xlabel="$\\mu_{E}$",fontsize=20)
ax_SEI_2D[0].set_ylabel(ylabel="$\\delta$",fontsize=20)
#ax_SEI_2D[0].set_title(label="Dynamics of outputs of SEI-model",fontsize=20)
ax_SEI_2D[0].xaxis.set_tick_params(labelsize=15)
ax_SEI_2D[0].yaxis.set_tick_params(labelsize=15)
#---------------------------------------------------------------------------------
# Subplot B: beta vs kI
#---------------------------------------------------------------------------------
for index in range(len(clr_vec)):
    #if index == 0 or index == 1:
        # Plot the parameters against eachother yeah?
    ax_SEI_2D[1].plot(p_hat[index][:,1],p_hat[index][:,8],color=clr_vec[index],linewidth=3.0, linestyle = "dashdot")        
# Set the y-limit
#ax_SEI[1].set_ylim([0, 2.5])
# Set a grid and define a legend
ax_SEI_2D[1].grid()
#ax_SEI_2D[1].legend(loc='best',prop={"size":20})
# Set the x-labels and y-labels
ax_SEI_2D[1].set_xlabel(xlabel="$\\beta$",fontsize=20)
ax_SEI_2D[1].set_ylabel(ylabel="$k_{I}$",fontsize=20)
#ax_SEI_2D[1].set_title(label="Dynamics of outputs of SEI-model",fontsize=20)
ax_SEI_2D[1].xaxis.set_tick_params(labelsize=15)
ax_SEI_2D[1].yaxis.set_tick_params(labelsize=15)
#---------------------------------------------------------------------------------
# Subplot C: epsilon vs delta
#---------------------------------------------------------------------------------
for index in range(len(clr_vec)):
    #if index == 0 or index == 1:
        # Plot the parameters against eachother yeah?
    ax_SEI_2D[2].plot(p_hat[index][:,2],p_hat[index][:,3],color=clr_vec[index],linewidth=3.0, linestyle = "dashdot")        
# Set the y-limit
#ax_SEI[1].set_ylim([0, 2.5])
# Set a grid and define a legend
ax_SEI_2D[2].grid()
#ax_SEI_2D[1].legend(loc='best',prop={"size":20})
# Set the x-labels and y-labels
ax_SEI_2D[2].set_xlabel(xlabel="$\\delta$",fontsize=20)
ax_SEI_2D[2].set_ylabel(ylabel="$\\epsilon$",fontsize=20)
#ax_SEI_2D[1].set_title(label="Dynamics of outputs of SEI-model",fontsize=20)
ax_SEI_2D[2].xaxis.set_tick_params(labelsize=15)
ax_SEI_2D[2].yaxis.set_tick_params(labelsize=15)
#---------------------------------------------------------------------------------
# Show plot and save figure
#---------------------------------------------------------------------------------    
plt.show()
# Title and saving the figure
f_SEI_2D.savefig('../Figures/SEI_parameter_plot_2D.png')




