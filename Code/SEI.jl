# SEI model with two outputs
using Latexify
using LaTeXStrings
using StructuralIdentifiability
# ============================================================================
# ODE model
# ============================================================================
ode_SEI = @ODEmodel(
    S'(t) = c - beta*S(t)*I(t) - mu_S*S(t),
    E'(t) = (1-eps)*beta*S(t)*I(t) - delta*E(t) -mu_E*E(t),
    I'(t) = eps*beta*S(t)*I(t) + delta*E(t) - mu_I*I(t),
    y1(t) = k_E*E(t),
    y2(t) = k_I*I(t)       
)
# Assess identifiability
println(assess_identifiability(ode_SEI))
# Print the equations 
println(first(values(find_ioequations(ode_SEI))))
# Print to a file
s = latexify(string(first(values(find_ioequations(ode_SEI)))))
open("../tex_files_with_output_ODEs/SEI.tex", "w") do file
    write(file, s)
end

