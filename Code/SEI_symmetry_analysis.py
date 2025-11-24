# =========================================================
# SEI symmetry analysis
# Date: 2025-11-24,
# Written by: Johannes
# Description: We solve the linearised symmetry conditions for the SEI model in order to find the parameter symmetries.
# =========================================================
# LOAD LIBRARIES
from sympy import * # For symbolic equations
# =========================================================
def main():
    #----------------------------------------------------------------------
    #----------------------------------------------------------------------
    # PART 1: DEFINING THE ODEs
    #----------------------------------------------------------------------    #---------------------------------------------------------------------
    # Open a file in which we save all equations in LaTeX
    f = open("./SEI_symmetry_analysis_document/notes.tex","w")
    # Define all parameters as symbols\n",
    c, beta, mu_S, mu_E, mu_I, delta, upsilon, kE, kI = symbols('c, beta, mu_S, mu_E, mu_I, delta, upsilon, kE, kI')
    # Define our independent variable time
    t = symbols('t')
    # Define our three dependent variables being the states
    S = Function('S')(t)
    E = Function('E')(t)
    I = Function('I')(t)
    # Now, we define our three ODEs
    # ODE for S
    ODE_S = c - beta*S*I - mu_S*S
    ODE_S = Eq(Derivative(S,t,1),ODE_S)
    f.write("ODE for $S$:")
    f.write(latex(ODE_S,mode='equation').replace("\\\\begin{equation}","\\\\begin{equation}\\n").replace("\\\\end{equation}","\\\\quad\\\\,,\\\\label{eq:ODE_S}\\n\\\\end{equation}\\n"))
    # ODE for E
    ODE_E = (1-upsilon)*beta*S*I - delta*E - mu_E*E
    ODE_E = Eq(Derivative(E,t,1),ODE_E)
    f.write("ODE for $E$:")
    f.write(latex(ODE_E,mode='equation').replace("\\\\begin{equation}","\\\\begin{equation}\\n").replace("\\\\end{equation}","\\\\quad\\\\,,\\\\label{eq:ODE_E}\\n\\\\end{equation}\\n"))
    # ODE for I
    ODE_I = upsilon*beta*S*I + delta*E - mu_I*I
    ODE_I = Eq(Derivative(I,t,1),ODE_I)
    f.write("ODE for $I$:")
    f.write(latex(ODE_I,mode='equation').replace("\\\\begin{equation}","\\\\begin{equation}\\n").replace("\\\\end{equation}","\\\\quad\\\\,.\\\\label{eq:ODE_I}\\n\\\\end{equation}\\n"))
    #----------------------------------------------------------------------
    #----------------------------------------------------------------------
    # PART 2: SETTING UP THE LINEARISED SYMMETRY CONDITIONS
    #----------------------------------------------------------------------    #----------------------------------------------------------------------    # Define major function for eta_S. This is the function we want to find yeah?
    eta_S = Function('eta_S')(t,S,E,I,c,beta,mu_S,mu_E,mu_I,delta,upsilon, kE,kI)
    

    

# RUN THE MAIN FUNCTION
if __name__ == "__main__":
    main()   
