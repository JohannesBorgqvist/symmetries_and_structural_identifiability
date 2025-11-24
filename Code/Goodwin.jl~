# Linear model with various outputs

using Latexify
using LaTeXStrings
using StructuralIdentifiability


# ============================================================================
# Case 1: y=k*u1
# ============================================================================
ode_case_1 = @ODEmodel(
    u1'(t) = a - b*u1(t),
    u2'(t) = b*u1(t)-c*u2(t),
    u3'(t) = c*u2(t),
    y(t) = k*u1(t)
)


println(assess_identifiability(ode_case_1))

# Not everything is identifiable, so we may wonder what are the identifiable functions

#println(find_identifiable_functions(ode_case_1, with_states = true))

println(first(values(find_ioequations(ode_case_1))))
# Print to a file
s = latexify(string(first(values(find_ioequations(ode_case_1)))))
open("../tex_files_with_output_ODEs/linear_case_1.tex", "w") do file
    write(file, s)
end


# ============================================================================
# Case 2: y=k*u2
# ============================================================================
ode_case_2 = @ODEmodel(
    u1'(t) = a - b*u1(t),
    u2'(t) = b*u1(t)-c*u2(t),
    u3'(t) = c*u2(t),
    y(t) = k*u2(t)
)


println(assess_identifiability(ode_case_2))


println(first(values(find_ioequations(ode_case_2))))
# Print to a file
s = latexify(string(first(values(find_ioequations(ode_case_2)))))
open("../tex_files_with_output_ODEs/linear_case_2.tex", "w") do file
    write(file, s)
end



# ============================================================================
# Case 3: y=k*u3
# ============================================================================
ode_case_3 = @ODEmodel(
    u1'(t) = a - b*u1(t),
    u2'(t) = b*u1(t)-c*u2(t),
    u3'(t) = c*u2(t),
    y(t) = k*u3(t)
)


println(assess_identifiability(ode_case_3))


println(first(values(find_ioequations(ode_case_3))))
# Print to a file
s = latexify(string(first(values(find_ioequations(ode_case_3)))))
open("../tex_files_with_output_ODEs/linear_case_3.tex", "w") do file
    write(file, s)
end



# ============================================================================
# Case 4: y=k*u1*u2
# ============================================================================
ode_case_4 = @ODEmodel(
    u1'(t) = a - b*u1(t),
    u2'(t) = b*u1(t)-c*u2(t),
    u3'(t) = c*u2(t),
    y(t) = k*u1(t)*u2(t)
)


println(assess_identifiability(ode_case_4))


println(first(values(find_ioequations(ode_case_4))))
# Print to a file
s = latexify(string(first(values(find_ioequations(ode_case_4)))))
open("../tex_files_with_output_ODEs/linear_case_4.tex", "w") do file
    write(file, s)
end
