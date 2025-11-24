# Linear model with various outputs

using Latexify
using LaTeXStrings
using StructuralIdentifiability


# ============================================================================
# Case 1: y=k*u1
# ============================================================================
ode_Goodwin = @ODEmodel(
    x1'(t) = -b * x1(t) + 1 / (c + x4(t)),
    x2'(t) = alpha * x1(t) - beta * x2(t),
    x3'(t) = gamma * x2(t) - delta * x3(t),
    x4'(t) = sigma * x4(t) * (gamma * x2(t) - delta * x3(t)) / x3(t),
    y(t) = k*x1(t)
)


println(assess_identifiability(ode_Goodwin))

# Not everything is identifiable, so we may wonder what are the identifiable functions

println(find_identifiable_functions(ode_Goodwin))
# Print to a file
s = latexify(string(first(values(find_identifiable_functions(ode_Goodwin)))))
open("../tex_files_with_output_ODEs/Goodwin_identifiable.tex", "w") do file
    write(file, s)
end


println(first(values(find_ioequations(ode_Goodwin))))
# Print to a file
#s = latexify(string(first(values(find_ioequations(ode_Goodwin)))))
s = string(first(values(find_ioequations(ode_Goodwin))))
open("../tex_files_with_output_ODEs/Goodwin.tex", "w") do file
    write(file, s)
end


