# Non-autonomous glucose model with a time-dependent output
using Latexify
using LaTeXStrings
using StructuralIdentifiability
# ============================================================================
# ODE model
# ============================================================================
ode_case_glucose = @ODEmodel(
    q1'(t) = u(t) +p1*q1(t)-p2*q2(t),
    q2'(t) = p3*q2(t)+p4*q1(t),
    y(t)=q1(t)/V
)
# Assess identifiability
println(assess_identifiability(ode_case_glucose))
# Print the equations 
println(first(values(find_ioequations(ode_case_glucose))))
# Print to a file
s = latexify(string(first(values(find_ioequations(ode_case_glucose)))))
open("../tex_files_with_output_ODEs/glucose.tex", "w") do file
    write(file, s)
end
