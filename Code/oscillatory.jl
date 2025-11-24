# Simple but tricky and overparametrised oscillatory system

using Latexify
using LaTeXStrings
using StructuralIdentifiability


# ============================================================================
# Oscillatory two state system
# ============================================================================
ode_oscillatory = @ODEmodel(
    u'(t) = a*b*v(t)-c,
    v'(t) = d*u(t)-c*d,
    y(t) = v(t)/e
)


println(assess_identifiability(ode_oscillatory))

# Not everything is identifiable, so we may wonder what are the identifiable functions

println(find_identifiable_functions(ode_oscillatory))
# Print to a file
s = latexify(string(first(values(find_identifiable_functions(ode_oscillatory)))))
open("../tex_files_with_output_ODEs/oscillatory_identifiable.tex", "w") do file
    write(file, s)
end


println(first(values(find_ioequations(ode_oscillatory))))
# Print to a file
#s = latexify(string(first(values(find_ioequations(ode_oscillatory)))))
s = string(first(values(find_ioequations(ode_oscillatory))))
open("../tex_files_with_output_ODEs/oscillatory.tex", "w") do file
    write(file, s)
end


