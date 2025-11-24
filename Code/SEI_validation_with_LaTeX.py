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
    # Unpack arbitrary coefficients
    alpha_1 = alpha[0]
    alpha_2 = alpha[1]
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
#==============================================================
# Function 3: plot_LaTeX_2D
# This function enables us to reproduce our plots using pgfplots in LaTeX in 2D
def plot_LaTeX_2D(t,y,file_str,plot_str,legend_str):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if len(legend_str)==0:
        temp_str = "\\addplot[\nforget plot,\n" + plot_str+ "\n]\n"
    else:
        temp_str = "\\addplot[\n" + plot_str+ "\n]\n"
    # Add the coordinates
    temp_str += "coordinates {%\n"
    # Loop over the input files and add
    # them to the file
    if len(t)==0:
        temp_str += "(" + str(t[0]) + "," + str(y[0]) + ")\n"
        # The plotting is done, let's close the shop    
        temp_str += "};\n"        
    else:
        for i in range(len(t)):
            temp_str += "(" + str(t[i]) + "," + str(y[i]) + ")\n"
        # The plotting is done, let's close the shop    
        temp_str += "};\n"
    # Add a legend if one is provided
    if len(legend_str) > 0:
        temp_str += "\\addlegendentry{" + legend_str + "}\n"
    # Finally, we write the huge string
    # we have created
    f.write("%s"%(temp_str))
    # Close the file
    f.close()
#==============================================================
# Function 4: plot_LaTeX_3D
# This function enables us to reproduce our plots using pgfplots in LaTeX in 3D
def plot_LaTeX_3D(x,y,z,file_str,plot_str,legend_str,surfaceNotCurve):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if surfaceNotCurve:
        if len(legend_str)==0:
            temp_str = "\\addplot3[forget plot," + plot_str+ "]\n"
        else:
            temp_str = "\\addplot3[" + plot_str+ "]\n"
    else:
        if len(legend_str)==0:
            temp_str = "\\addplot3+[forget plot," + plot_str+ "]\n"
        else:
            temp_str = "\\addplot3+[" + plot_str+ "]\n"        
    # Add the coordinates
    temp_str += "coordinates {%\n"
    if len(x)==0:
        temp_str += "(" + str(x[0]) + "," + str(y[0]) + "," + str(z[0])+ ")\n"
        # The plotting is done, let's close the shop    
        temp_str += "};\n"        
    else:
        for i in range(len(x)):
            temp_str += "(" + str(x[i]) + "," + str(y[i]) + "," + str(z[i])+ ")\n"
        # The plotting is done, let's close the shop    
        temp_str += "};\n"
    # # Loop over the input files and add
    # # them to the file
    # for index in range(len(data)):
    #     temp_str += "("+str(data[index][0]) + "," + str(data[index][1]) + "," + str(data[index][2]) + ")"
    #     if index>0:
    #         if index < len(data)-1 and data[index][1] < data[index+1][1]:
    #             temp_str += "\n"
    #         elif index == len(data)-1:
    #             temp_str += "\n"
    #         else:
    #             temp_str += "  "
    # # The plotting is done, let's close the shop    
    # temp_str += "};\n"
    # Add a legend if one is provided
    if len(legend_str) > 0:
        temp_str += "\\addlegendentry{" + legend_str + "}\n"
    # Finally, we write the huge string
    # we have created
    f.write("%s"%(temp_str))
    # Close the file
    f.close()    
# =========================================================
# Solve ODE system
# =========================================================
# State the number of data point
num_data_points = 200
# Simulate time forwards
t = np.linspace(0, 100, num_data_points)
# Define parameters
p = [200, 1.0, 0.1, 0.05, 0.005, 0.005, 0.01, 0.5, 0.5]
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
# =========================================================
# Transform parameters and solve system again
# =========================================================
# Define a var_eps vector
var_eps_vec = np.linspace(0,1,100)
# Allocate a vector for the alpha values
alpha = [2,1]
# # Allocate memory for our SSs
# SS_yE = []
# SS_yI = []
# for var_eps in var_eps_vec:
#     # Define a var_eps vector
#     var_eps_vec_temp = np.linspace(0,var_eps,1000000)
#     # Solve ODE system for the parameters
#     p_hat_sys = odeint(param_symmetry_SEI, p, var_eps_vec_temp, args=(alpha,))
#     # Update the transformed parameter vector
#     p_hat = p_hat_sys[-1]
#     # Define the ODEs
#     SEI_states = odeint(SEI_outputs, y0, t, args=(p_hat,))
#     # Extract states
#     yI_hat, yE_hat, dot_yI_hat = SEI_states.T
#     # Calculate residuals
#     SS_yE.append(np.inner(yE_hat-yE,yE_hat-yE)/len(yI))
#     SS_yI.append(np.inner(yI_hat-yI,yI_hat-yI)/len(yI))
# # Re-cast lists for the outputs yI and yE as arrays
# SS_yE = np.array(SS_yE)
# SS_yI = np.array(SS_yI)
# #---------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------
# # Plot SEI solutions
# #---------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------
# #Define the first figure
# f_SEI, ax_SEI = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
# #---------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------
# # Subplot A
# #---------------------------------------------------------------------------------
# # Plot the simulated data and the underlying model
# ax_SEI[0].plot(t,yE,color=(0/256,68/256,27/256),label="$y_{E}(t)$",linewidth=3.0)
# ax_SEI[0].plot(t,yI,color=(35/256,139/256,69/256),label="$y_{I}(t)$",linewidth=3.0)
# # Set the y-limit
# #ax_SEI[0].set_ylim([0, 2.5])
# # Set a grid and define a legend
# ax_SEI[0].grid()
# ax_SEI[0].legend(loc='best',prop={"size":20})
# # Set the x-labels and y-labels
# ax_SEI[0].set_xlabel(xlabel="Time, $t$",fontsize=20)
# ax_SEI[0].set_ylabel(ylabel="Population density",fontsize=20)
# ax_SEI[0].set_title(label="Dynamics of outputs of SEI-model",fontsize=20)
# ax_SEI[0].xaxis.set_tick_params(labelsize=15)
# ax_SEI[0].yaxis.set_tick_params(labelsize=15)
# #---------------------------------------------------------------------------------
# # Subplot B
# #---------------------------------------------------------------------------------
# # Plot the simulated data and the underlying model
# ax_SEI[1].plot(var_eps_vec,SS_yE,color=(2/256,56/256,88/256),label="$\\sqrt{\\frac{\\sum_{i=1}^{n}\\left(y_{E}(t_{i},\\vec{\\theta})-y_{E}(t_{i},\\hat{\\vec{\\theta}}(\\varepsilon))\\right)^{2}}{n}}$",linewidth=3.0)
# ax_SEI[1].plot(var_eps_vec,SS_yI,color=(35/256,139/256,69/256),label="$\\sqrt{\\frac{\\sum_{i=1}^{n}\\left(y_{I}(t_{i},\\vec{\\theta})-y_{I}(t_{i},\\hat{\\vec{\\theta}}(\\varepsilon))\\right)^{2}}{n}}$",linewidth=3.0)
# # # Set a grid and define a legend
# ax_SEI[1].grid()
# ax_SEI[1].legend(loc='best',prop={"size":20})
# # Set the x-labels and y-labels
# ax_SEI[1].set_xlabel(xlabel="$\\varepsilon$",fontsize=20)
# ax_SEI[1].set_ylabel(ylabel="Sum of squares, $\\mathrm{SS}(\\varepsilon)$",fontsize=20)
# ax_SEI[1].set_title(label="The fit is invariant under transformation by $\\Gamma_{\\varepsilon}^{\\vec{\\theta}}$",fontsize=20)
# ax_SEI[1].xaxis.set_tick_params(labelsize=15)
# ax_SEI[1].yaxis.set_tick_params(labelsize=15)
# plt.show()
# # Title and saving the figure
# f_SEI.savefig('../Figures/SEI_numerical_validation.png')


# # Figure path
# #figure_path = "/home/johborgq/Dropbox/Work/fancy_research_images/SEI_numerical_validation_of_parameter_symmetry"
# figure_path = "/home/johannes/Dropbox/Work/fancy_research_images/SEI_numerical_validation_of_parameter_symmetry"
# # Save things in LaTeX as well
# plot_LaTeX_2D(t,yE,figure_path + "/Input/outputs.tex","color=clr_1,line width=1.5pt,","$y_{\mathrm{E}}(t)$")
# plot_LaTeX_2D(t,yI,figure_path + "/Input/outputs.tex","color=clr_2,line width=1.5pt,","$y_{\mathrm{I}}(t)$")
# #plot_LaTeX_2D(var_eps_vec,SS_yE,figure_path + "/Input/SS.tex","color=clr_1,line width=1.5pt,","$\\sqrt{\\frac{\\sum_{i=1}^{n}\\left(y_{E}(t_{i},\\vec{\\theta})-y_{E}(t_{i},\\hat{\\vec{\\theta}}(\\varepsilon))\\right)^{2}}{n}}$")
# plot_LaTeX_2D(var_eps_vec,SS_yE,figure_path + "/Input/SS.tex","color=clr_1,line width=1.5pt,","$\\frac{1}{T}\\bigintsss_{0}^{T}\\left(y_{E}(t,\\vec{\\theta})-y_{E}(t,\\hat{\\vec{\\theta}}(\\varepsilon))\\right)^{2}\\mathrm{d}t$")
# #plot_LaTeX_2D(var_eps_vec,SS_yI,figure_path + "/Input/SS.tex","color=clr_2,line width=1.5pt,","$\\sqrt{\\frac{\\sum_{i=1}^{n}\\left(y_{I}(t_{i},\\vec{\\theta})-y_{I}(t_{i},\\hat{\\vec{\\theta}}(\\varepsilon))\\right)^{2}}{n}}$")
# plot_LaTeX_2D(var_eps_vec,SS_yI,figure_path + "/Input/SS.tex","color=clr_2,line width=1.5pt,","$\\frac{1}{T}\\bigintsss_{0}^{T}\\left(y_{I}(t,\\vec{\\theta})-y_{I}(t,\\hat{\\vec{\\theta}}(\\varepsilon))\\right)^{2}\\mathrm{d}t$")



# =========================================================
# Illustrate parameter symmetries
# =========================================================
# Define a new slightly longer var_eps vector
#var_eps_vec = np.linspace(0,5,500) # For the 2D plots
var_eps_vec = np.linspace(0,3,500) # For the 3D plots
# Allocate a parameter vector for the symmetries
p_vec = [p]
# Fix the seed of the generator so we get the same result each time
np.random.seed(55)
# Add five more parameter vectors, shall we?
for i in range(5):
    p_vec.append(np.array(np.multiply(np.random.rand(1,len(p)),p)[0]))
# Ok, let's use our var_eps vec to define some transformed parameters shall we?
p_trans = []
# Loop over the original parameters and calculate the transformed parameter vectors
for p_temp in p_vec:
     p_trans.append(odeint(param_symmetry_SEI, p_temp, var_eps_vec, args=(alpha,)))
# Plot a single point, yeah?
# Define a figure path
figure_path = "/home/johannes/Dropbox/Work/fancy_research_images/SEI_parameter_symmetries_2D"
# ------------------------------------------------------------------------------------------------------------------------------------------------
# Parameter pairs in 2D
# ------------------------------------------------------------------------------------------------------------------------------------------------
# Define all indices of pairwise coupled parameters
pair_indices = [(5,2), (8,1), (2,3)]
# Define all file names
file_names = ["muE_vs_delta", "kI_vs_beta", "delta_vs_epsilon"]
# Loop over file names and pair indices and save in lovely files
for index, file_name in enumerate(file_names):
    # Extract indices yeah?
    i, j = pair_indices[index]
    # Loop over number of random numbers, yeah?
    for inner_index in range(len(p_vec)):
        # Start with delta vs muE
        x = p_trans[inner_index][:,i]
        y = p_trans[inner_index][:,j]
        # Save the figures
        plot_LaTeX_2D([x[0]],[y[0]],figure_path + "/Input/" + file_name + ".tex","color=clr_" + str(inner_index+1) + ",line width=7pt,mark=diamond*",[])
        # Define a previous index
        prev_index = 0
        number_of_points = 10
        number_of_iter = int(len(var_eps_vec)/number_of_points)
        # Loop to plot an arrow
        #for i in range(number of points):
        for repeating in range(number_of_iter):
            if prev_index == 0:
                plot_LaTeX_2D(x[prev_index:prev_index+number_of_points],y[prev_index:prev_index+number_of_points],figure_path + "/Input/" + file_name + ".tex","color=clr_" + str(inner_index+1) + ",line width=1.5pt,thick,dashed,-latex,->",[])
            else:
                plot_LaTeX_2D(x[prev_index-1:prev_index-1+number_of_points],y[prev_index-1:prev_index-1+number_of_points],figure_path + "/Input/" + file_name + ".tex","color=clr_" + str(inner_index+1) + ",line width=1.5pt,thick,dashed,-latex,->",[])
            prev_index += number_of_points
    
# ------------------------------------------------------------------------------------------------------------------------------------------------
# Parameter triplets in 3D
# ------------------------------------------------------------------------------------------------------------------------------------------------
# Define a figure path
figure_path = "/home/johannes/Dropbox/Work/fancy_research_images/SEI_parameter_symmetries_3D"
# Define all indices of pairwise coupled parameters
triplet_indices = [(1,3,0), (1,2,7)]
# Define all file names
file_names = ["beta_vs_epsilon_vs_c_3D", "beta_vs_delta_vs_kE_3D"]
file_names_2D = ["beta_vs_epsilon_vs_c_2D", "beta_vs_delta_vs_kE_2D"]
# Loop over file names and pair indices and save in lovely files
for index, file_name in enumerate(file_names):
    # Extract indices yeah?
    i, j, k = triplet_indices[index]
    # Loop over number of random numbers, yeah?
    for inner_index in range(len(p_vec)):
        # Extract the data we want to plot yeah?
        x = [p_trans[inner_index][0,i]]
        y = [p_trans[inner_index][0,j]]
        z = [p_trans[inner_index][0,k]]
        # Save the figures
        plot_LaTeX_3D(x,y,z,figure_path + "/Input/" + file_name + ".tex","color=clr_" + str(inner_index+1) + ",line width=7pt,mark=diamond*",[],False)
        # Plot the 2D version as well
        if index == 0:
            plot_LaTeX_2D(x,[y[0]*z[0]],figure_path + "/Input/" + file_names_2D[index] + ".tex","color=clr_" + str(inner_index+1) + ",line width=7pt,mark=diamond*",[])
        else:
            plot_LaTeX_2D(x,[y[0]/z[0]],figure_path + "/Input/" + file_names_2D[index] + ".tex","color=clr_" + str(inner_index+1) + ",line width=7pt,mark=diamond*",[])        
        # Define a previous index
        prev_index = 0
        number_of_points = 10
        number_of_iter = int(len(var_eps_vec)/number_of_points)
        # Loop to plot an arrow
        #for i in range(number of points):
        for repeating in range(number_of_iter):
            # Extract data points yeah?
            x = p_trans[inner_index][:,i]
            y = p_trans[inner_index][:,j]
            z = p_trans[inner_index][:,k]            
            if prev_index == 0:
                plot_LaTeX_3D(x[prev_index:prev_index+number_of_points],y[prev_index:prev_index+number_of_points],z[prev_index:prev_index+number_of_points],figure_path + "/Input/" + file_name + ".tex","color=clr_" + str(inner_index+1) + ",line width=1.5pt,thick,dashed,-latex,->,mark=none",[],False)
                if index == 0:
                    plot_LaTeX_2D(x[prev_index:prev_index+number_of_points],y[prev_index:prev_index+number_of_points]*z[prev_index:prev_index+number_of_points],figure_path + "/Input/" + file_names_2D[index] + ".tex","color=clr_" + str(inner_index+1) + ",line width=1.5pt,thick,dashed,-latex,->,mark=none",[])
                else:
                    plot_LaTeX_2D(x[prev_index:prev_index+number_of_points],y[prev_index:prev_index+number_of_points]//z[prev_index:prev_index+number_of_points],figure_path + "/Input/" + file_names_2D[index] + ".tex","color=clr_" + str(inner_index+1) + ",line width=1.5pt,thick,dashed,-latex,->,mark=none",[])                        
            else:
                plot_LaTeX_3D(x[prev_index-1:prev_index-1+number_of_points],y[prev_index-1:prev_index-1+number_of_points],z[prev_index-1:prev_index-1+number_of_points],figure_path + "/Input/" + file_name + ".tex","color=clr_" + str(inner_index+1) + ",line width=1.5pt,thick,dashed,-latex,->,mark=none",[],False)
                if index == 0:
                    plot_LaTeX_2D(x[prev_index-1:prev_index-1+number_of_points],y[prev_index-1:prev_index-1+number_of_points]*z[prev_index-1:prev_index-1+number_of_points],figure_path + "/Input/" + file_names_2D[index] + ".tex","color=clr_" + str(inner_index+1) + ",line width=1.5pt,thick,dashed,-latex,->,mark=none",[])
                else:
                    plot_LaTeX_2D(x[prev_index-1:prev_index-1+number_of_points],np.divide(y[prev_index-1:prev_index-1+number_of_points],z[prev_index-1:prev_index-1+number_of_points]),figure_path + "/Input/" + file_names_2D[index] + ".tex","color=clr_" + str(inner_index+1) + ",line width=1.5pt,thick,dashed,-latex,->,mark=none",[])                        
            prev_index += number_of_points
