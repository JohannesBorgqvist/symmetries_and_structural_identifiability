# =========================================================
# SEI symmetry analysis
# Date: 2025-11-24,
# Written by: Johannes
# Description: We solve the linearised symmetry conditions for the SEI model in order to find the parameter symmetries.
# =========================================================
# LOAD LIBRARIES
from sympy import * # For symbolic equations
# For manipulating monomials
from sympy.polys.monomials import itermonomials
from sympy.polys.orderings import monomial_key
# =========================================================
def main():
    #----------------------------------------------------------------------
    #----------------------------------------------------------------------
    # PART 1: DEFINING THE ODEs
    #----------------------------------------------------------------------
    #---------------------------------------------------------------------
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
    #----------------------------------------------------------------------
    #----------------------------------------------------------------------
    # Define major function for eta_S. This is the function we want to find yeah?
    eta_S = Function('eta_S')(t,S,E,I,c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    # Define parameter infinitesimals
    chi_c=Function('chi_c')(c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    chi_beta=Function('chi_beta')(c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    chi_muS=Function('chi_muS')(c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    chi_muE=Function('chi_muE')(c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    chi_muI=Function('chi_muI')(c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    chi_delta=Function('chi_delta')(c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    chi_upsilon=Function('chi_upsilon')(c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    chi_kE=Function('chi_kE')(c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    chi_kI=Function('chi_kI')(c,beta,mu_S,mu_E,mu_I,delta,upsilon,kE,kI)
    # Define a vector with the infinitesimals
    chi_vec = Matrix(9,1,[chi_c, chi_beta, chi_muS, chi_muE, chi_muI, chi_delta, chi_upsilon, chi_kE, chi_kI])
    # Define our new infinitesimals
    eta_E = -(chi_kE/kE)*E
    eta_I = -(chi_kI/kI)*I
    # Print this to our lovely notes
    f.write("Equation for $\\eta_{E}$:")
    f.write("\\begin{equation}\n\\eta_{E}=%s\\,\\label{eq:eta_E_EQ}\n\\end{equation}"%(latex(eta_E)))
    f.write("Equation for $\\eta_{I}$:")
    f.write("\\begin{equation}\n\\eta_{I}=%s\\,\\label{eq:eta_I_EQ}\n\\end{equation}"%(latex(eta_I)))    
    #----------------------------------------------------------------------------------------------------
    # Let's set up the linearised symmetry condition for E
    lin_sym_E_lhs = Derivative(eta_E,t,1).doit().subs(ODE_E.lhs,ODE_E.rhs)
    # Let's put up the right hand side in parts please
    lin_sym_E_rhs = Derivative(ODE_E.rhs,S).doit()*eta_S+Derivative(ODE_E.rhs,E).doit()*eta_E+Derivative(ODE_E.rhs,I).doit()*eta_I
    lin_sym_E_rhs += Derivative(ODE_E.rhs,c).doit()*chi_c+Derivative(ODE_E.rhs,beta).doit()*chi_beta+Derivative(ODE_E.rhs,delta).doit()*chi_delta
    lin_sym_E_rhs += Derivative(ODE_E.rhs,upsilon).doit()*chi_upsilon+Derivative(ODE_E.rhs,kE).doit()*chi_kE+Derivative(ODE_E.rhs,kI).doit()*chi_kI
    lin_sym_E_rhs += Derivative(ODE_E.rhs,mu_S).doit()*chi_muS+Derivative(ODE_E.rhs,mu_E).doit()*chi_muE+Derivative(ODE_E.rhs,mu_I).doit()*chi_muI
    # Add the lhs and the rhs to get the lin_sym for E
    lin_sym_E = Eq(lin_sym_E_lhs,lin_sym_E_rhs)
    f.write("Linearised symmetry condition for $E$:")
    f.write(latex(lin_sym_E,mode='equation').replace("\\\\begin{equation}","\\\\begin{equation}\\n").replace("\\\\end{equation}","\\\\quad\\\\,.\\\\label{eq:lin_sym_E_sympy}\\n\\\\end{equation}\\n"))
    #----------------------------------------------------------------------------------------------------
    # Let's set up the linearised symmetry condition for I
    lin_sym_I_lhs = Derivative(eta_I,t,1).doit().subs(ODE_I.lhs,ODE_I.rhs)
    # Let's put up the right hand side in parts please
    lin_sym_I_rhs = Derivative(ODE_I.rhs,S).doit()*eta_S+Derivative(ODE_I.rhs,E).doit()*eta_E+Derivative(ODE_I.rhs,I).doit()*eta_I
    lin_sym_I_rhs += Derivative(ODE_I.rhs,c).doit()*chi_c+Derivative(ODE_I.rhs,beta).doit()*chi_beta+Derivative(ODE_I.rhs,delta).doit()*chi_delta
    lin_sym_I_rhs += Derivative(ODE_I.rhs,upsilon).doit()*chi_upsilon+Derivative(ODE_I.rhs,kE).doit()*chi_kE+Derivative(ODE_I.rhs,kI).doit()*chi_kI
    lin_sym_I_rhs += Derivative(ODE_I.rhs,mu_S).doit()*chi_muS+Derivative(ODE_I.rhs,mu_E).doit()*chi_muE+Derivative(ODE_I.rhs,mu_I).doit()*chi_muI
    # Add the lhs and the rhs to get the lin_sym for I
    lin_sym_I = Eq(lin_sym_I_lhs,lin_sym_I_rhs)
    f.write("Linearised symmetry condition for $I$:")
    f.write(latex(lin_sym_I,mode='equation').replace("\\\\begin{equation}","\\\\begin{equation}\\n").replace("\\\\end{equation}","\\\\quad\\\\,.\\\\label{eq:lin_sym_I_sympy}\\n\\\\end{equation}\\n"))
    #----------------------------------------------------------------------------------------------------
    # Solve the linearised symmetry condition for E for the factor "(1-upsilon)*beta*I*eta_S"
    eq_lin_sym_E = Eq((1-upsilon)*beta*I*eta_S,expand(simplify(solve(lin_sym_E.lhs-lin_sym_E.rhs,(1-upsilon)*beta*I*eta_S)[0])))
    f.write("\n\nEquation for factor ``$(1-\\upsilon)\\beta{I}\\eta_{S}$'' from lin.sym for $E$:\n")
    f.write(latex(eq_lin_sym_E,mode='equation').replace("\\\\begin{equation}","\\\\begin{equation}\\n").replace("\\\\end{equation}","\\\\quad\\\\,.\\\\label{eq:lin_sym_E_factor_sympy}\\n\\\\end{equation}\\n"))

    #----------------------------------------------------------------------------------------------------
    # Solve both linearised symmetry condition for eta_S
    eta_S_1 = solve(lin_sym_E.lhs-lin_sym_E.rhs,eta_S)[0]
    eta_S_2 = solve(lin_sym_I.lhs-lin_sym_I.rhs,eta_S)[0]
    # Equate these conditions and simplify
    eta_S_equality = simplify(eta_S_1 - eta_S_2)
    eta_S_equality, d = fraction(eta_S_equality)
    eta_S_equality = Eq(expand(eta_S_equality),0)
    # Save this equality in our notes
    f.write("\n\nEquality stemming from $\\eta$_{S}:\n\n")
    f.write(latex(eta_S_equality,mode='equation'))
    #----------------------------------------------------------------------------------------------------
    # Decompose the linearised symmetry condition with respect to S using iterative monomials of the states S, E and I.
    #----------------------------------------------------------------------------------------------------
    # States
    states = [S,E,I]
    # Create iterated monomials using these states
    monomials = sorted(itermonomials(states, 10),key=monomial_key('grlex', states))
    # Allocate our determining equations
    det_eqs = []
    # Allocate a temporary equation
    temp_eq = 0
    # Loop over monomials and extract determining equations
    for monomial in monomials:
        # Extract a symmetry condition
        temp_eq = eta_S_equality.lhs
        # Special case when the monomial is 1
        if monomial==1:
            # For this special case we save only the constant term
            # by setting all derivatives and states to zero
            for state in states:
                temp_eq = temp_eq.subs(state,0).doit()
        else:
            # Extract the coefficient
            temp_eq = temp_eq.coeff(monomial)
            # Just to be sure we only save the constant thingy
            for state in states:
                temp_eq = temp_eq.subs(state,0).doit()
        # Save all non-zero determining equations
        if temp_eq != 0:
            det_eqs.append((monomial,temp_eq))    
    #----------------------------------------------------------------------------------------------------
    # Save these equations in our notes
    #----------------------------------------------------------------------------------------------------
    # Add an explanatory text for the text string for the determining equation     
    latex_str = "\n\nDetermining equations stemming from the linearised symmetry condition for S:\n\n"
    # We aling equations
    latex_str += "\\begin{align}\n"
    # Loop over the determining equations and add them to our string
    # Loop over the determining equations and print them to our LaTeX string
    for index,det_eq in enumerate(det_eqs):
        latex_str += latex(det_eq[0]) + ":&\\quad" + latex(det_eq[1]) + "&=0\\,,\\label{eq:det_eq_" + str(index) + "}\\\\\n"    # End the aligning of equations
    latex_str += "\n\\end{align}\n"
    # Write it to our file
    f.write(latex_str)    
    #----------------------------------------------------------------------------------------------------
    # Setup our coefficient matrix
    #----------------------------------------------------------------------------------------------------
    # We begin by assembling a small vector of parameter infinitesimals
    chi_0 =[chi_upsilon, chi_kI, chi_kE, chi_muE, chi_muI, chi_delta]
    # Assemble matrix
    matrix_list = []
    # Loop over the determining equations and save them in our matrix list
    for det_eq in det_eqs:
        for chi_val in chi_0:
            matrix_list.append(det_eq[1].coeff(chi_val))
    # Assemble the matrix using the matrix list
    mat_chi = Matrix(len(det_eqs),len(chi_0),matrix_list)
    # Write this to our notes
    f.write("\n\nCoefficient matrix:")
    f.write("$$" + latex(mat_chi) + "\\,.$$\n\n")
    #----------------------------------------------------------------------------------------------------
    # Calculate the nullspace of our beloved matrix
    #----------------------------------------------------------------------------------------------------
    chi_nullspace = mat_chi.nullspace()
    f.write("\n\nNull space of matrix:\n")
    f.write("$$\\mathcal{N}(M)=" + latex(chi_nullspace) + "\\,.$$\n")
    #----------------------------------------------------------------------------------------------------
    # Calculate the solution
    #----------------------------------------------------------------------------------------------------
    # Allocate three arbitrary parameters
    a1, a2, a3 = symbols('a1 a2 a3')
    # Express solution as a linear combination of the basis vectors of the null space
    chi_0_sol = Matrix(6,1,list(simplify(a1*kE*chi_nullspace[0]+a2*delta*chi_nullspace[1]+a3*delta*chi_nullspace[2])))
    # Make this into an equation which we can print to a file
    f.write("\n\nSolution to matrix system:\n")
    f.write(latex(chi_0_sol,mode='equation'))
    #----------------------------------------------------------------------------------------------------
    # Combine these last equations in order to get an equation for eta_S
    #----------------------------------------------------------------------------------------------------    
    eta_S_particular = eta_S_1
    for i in range(len(chi_0_sol)):
        eta_S_particular = eta_S_particular.subs(chi_0[i],chi_0_sol[i])
    eta_S_particular = expand(simplify(eta_S_particular))
    f.write("\n\nEquation for factor $\\eta_{S}$:\n")
    f.write("$$\\eta_{S}="+latex(eta_S_particular)+"\\,.$$\\n")
    #----------------------------------------------------------------------------------------------------
    # Let's set up the linearised symmetry condition for S.
    # For the LHS we differentiate with respect to t and then we
    # substitute the RHS for the ODE. 
    lin_sym_S_lhs = Derivative(eta_S_particular,t,1).doit()
    lin_sym_S_lhs = lin_sym_S_lhs.subs(ODE_S.lhs,ODE_S.rhs)
    lin_sym_S_lhs = lin_sym_S_lhs.subs(ODE_E.lhs,ODE_E.rhs)
    lin_sym_S_lhs = lin_sym_S_lhs.subs(ODE_I.lhs,ODE_I.rhs)
    # Let's put up the right hand side in parts please
    lin_sym_S_rhs = Derivative(ODE_S.rhs,S).doit()*eta_S_particular+Derivative(ODE_S.rhs,E).doit()*eta_E+Derivative(ODE_S.rhs,I).doit()*eta_I
    lin_sym_S_rhs += Derivative(ODE_S.rhs,c).doit()*chi_c+Derivative(ODE_S.rhs,beta).doit()*chi_beta+Derivative(ODE_S.rhs,delta).doit()*chi_delta
    lin_sym_S_rhs += Derivative(ODE_S.rhs,upsilon).doit()*chi_upsilon+Derivative(ODE_S.rhs,kE).doit()*chi_kE+Derivative(ODE_S.rhs,kI).doit()*chi_kI
    lin_sym_S_rhs += Derivative(ODE_S.rhs,mu_S).doit()*chi_muS+Derivative(ODE_S.rhs,mu_E).doit()*chi_muE+Derivative(ODE_S.rhs,mu_I).doit()*chi_muI
    # Add the lhs and the rhs to get the lin_sym for I
    lin_sym_temp = lin_sym_S_lhs-lin_sym_S_rhs
    # Loop through the parameter infinitesimals and substitute their values
    for i in range(len(chi_0)):
        lin_sym_temp = lin_sym_temp.subs(chi_0[i],chi_0_sol[i])
    lin_sym_temp = simplify(lin_sym_temp)
    lin_sym_temp, d = fraction(lin_sym_temp)
    lin_sym_temp = expand(lin_sym_temp)
    lin_sym_S = Eq(lin_sym_temp,0)
    f.write("\\nLinearised symmetry condition for $S$:\\n\\n")
    f.write(latex(lin_sym_S,mode='equation').replace("\\\\begin{equation}","\\\\begin{equation}\\n").replace("\\\\end{equation}","\\\\quad\\\\,.\\\\label{eq:lin_sym_S_sympy}\\n\\\\end{equation}\\n"))
    #----------------------------------------------------------------------------------------------------
    # Decompose the linearised symmetry condition with respect to S using iterative monomials of the states S, E and I.
    # States
    states = [S,E,I]
    # Create iterated monomials using these states
    monomials = sorted(itermonomials(states, 10),key=monomial_key('grlex', states))
    # Allocate our determining equations
    det_eqs = []
    # Allocate a temporary equation
    temp_eq = 0
    # Loop over monomials and extract determining equations
    for monomial in monomials:
        # Extract a symmetry condition
        temp_eq = lin_sym_S.lhs
        # Special case when the monomial is 1
        if monomial==1:
            # For this special case we save only the constant term
            # by setting all derivatives and states to zero
            for state in states:
                temp_eq = temp_eq.subs(state,0).doit()
        else:
            # Extract the coefficient
            temp_eq = temp_eq.coeff(monomial)
            # Just to be sure we only save the constant thingy
            for state in states:
                temp_eq = temp_eq.subs(state,0).doit()
        # Save all non-zero determining equations
        if temp_eq != 0:
            det_eqs.append((monomial,temp_eq))      
    #----------------------------------------------------------------------------------------------------
    # Print all determining equations in a LaTeX string
    #----------------------------------------------------------------------------------------------------
    # Add an explanatory text for the text string for the determining equation     
    latex_str = "\n\nDetermining equations stemming from the linearised symmetry condition for S:\n\n"
    # We aling equations
    latex_str += "\\begin{align}\n"
    # Loop over the determining equations and add them to our string
    # Loop over the determining equations and print them to our LaTeX string
    for index,det_eq in enumerate(det_eqs):
        latex_str += latex(det_eq[0]) + ":&\\quad" + latex(det_eq[1]) + "&=0\\,,\\label{eq:det_eq_" + str(index) + "}\\\\\n"    # End the aligning of equations
    latex_str += "\n\\end{align}\n"
    # Write it to our file
    f.write(latex_str)
    #----------------------------------------------------------------------------------------------------
    # Update the determining equations with the substitution a3=-a2
    #----------------------------------------------------------------------------------------------------
    det_eqs_new = []
    temp_eq = 0
    for det_eq in det_eqs:
        temp_eq = simplify(det_eq[1].subs(a3,-a2))
        if temp_eq != 0:
            det_eqs_new.append((det_eq[0],temp_eq))
    #----------------------------------------------------------------------------------------------------
    # Print all these determining equations after this alph substitution a in a LaTeX string
    #----------------------------------------------------------------------------------------------------
    # Add an explanatory text for the text string for the determining equation     
    latex_str = "\n\nDetermining equations stemming from the linearised symmetry condition for S:\n\n"
    # We aling equations
    latex_str += "\\begin{align}\n"
    # Loop over the determining equations and print them to our LaTeX string
    for index,det_eq in enumerate(det_eqs_new):
        latex_str += latex(det_eq[0]) + ":&\\quad" + latex(det_eq[1]) + "&=0\\,,\\label{eq:det_eq_" + str(index) + "}\\\\\n"    # End the aligning of equations
    latex_str += "\n\\end{align}\n"
    # Write it to our file
    f.write(latex_str)
    #----------------------------------------------------------------------------------------------------
    # Do some final manipulations to find the two last infinitesimals. 
    #----------------------------------------------------------------------------------------------------
    # Extract the last equation for finding the parameter infinitesimal for beta
    det_eq_for_beta = Eq(chi_beta,simplify(solve(det_eqs_new[2][1],chi_beta)[0]))
    f.write("\n\nEquation for $\\chi_{\\beta}$:")
    f.write(latex(det_eq_for_beta,mode='equation'))
    # Extract the first equation for finding the parameter infinitesimal for c
    det_eq_for_c = Eq(chi_c,simplify(solve(det_eqs_new[0][1].subs(det_eq_for_beta.lhs,det_eq_for_beta.rhs),chi_c)[0]))
    f.write("\n\nEquation for $\\chi_{c}$:")
    f.write(latex(det_eq_for_c,mode='equation'))    
    #----------------------------------------------------------------------------------------------------
    # Finally, save the solutions in a chi_vec 
    #----------------------------------------------------------------------------------------------------
    # Allocate the final answer for the parameter infinitesimals
    chi_answer = Matrix(9,1,list(chi_vec))
    # Loop through chi_vec and substitue the zero infinitesimals
    for i in range(len(chi_answer)):
        if chi_answer[i] == chi_muI:
            chi_answer[i] = 0
        elif chi_answer[i] == chi_muS:
            chi_answer[i] = 0
        elif chi_answer[i] == det_eq_for_beta.lhs:
            chi_answer[i] = det_eq_for_beta.rhs
        elif chi_answer[i] == det_eq_for_c.lhs:
            chi_answer[i] = det_eq_for_c.rhs
    # Now, we do the final substitutions for chi_0
    for i in range(len(chi_answer)):
        for j in range(len(chi_0)):
            if chi_answer[i] == chi_0[j]:
                chi_answer[i] = simplify(chi_0_sol[j].subs(a3,-a2))
    # Write the final infinitesimals to a file
    f.write("\n\nThe parameter infinitesimals are:\n")
    f.write("\\begin{equation}\n" + latex(chi_vec) + "=" + latex(chi_answer) + "\\,.\\label{eq:parameter_infinitesimals}\n\\end{equation}\n\n")
    # We also print the updated version of the infinitesimal eta_S
    eta_S_final = Eq(eta_S,simplify(eta_S_particular.subs(a3,-a2).subs(det_eq_for_beta.lhs,det_eq_for_beta.rhs)))
    f.write("\n\nThe parameter infinitesimal $\\eta_{S}$:")
    f.write(latex(eta_S_final,mode='equation'))
    #----------------------------------------------------------------------------------------------------
    # Final prompt to the user 
    #----------------------------------------------------------------------------------------------------
    print("\n\n\t\tAll symmetry calculations for the SEI model were successful, and the results are saved in the \"notes.tex\" file. Good bye!\n\n ")
# RUN THE MAIN FUNCTION
if __name__ == "__main__":
    main()   
