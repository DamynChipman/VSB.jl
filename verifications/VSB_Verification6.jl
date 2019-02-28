# ===== VSB Verifications: Vortex Sheet Calculations =====

println("=== BEGINNING Verifications6.jl ===")

# === Imports ===
println("   IMPORTING PACKAGES...")
using VSB
using SimpleVPM
using LinearAlgebra
using Plots
pyplot()
using LaTeXStrings
#using Printf

cd("/Users/Damyn/Documents/BYU/FLOW Lab/VSB/verifications")

# === Boundary Geometry and Plotting Grids ===
println("   GENERATING GEOMETRY...")
filename = "naca0012.txt"

X_coor = []
Y_coor = []

open(filename, "r") do file
    line_number = 1
    for line in eachline(file)
        if line_number == 1 # Header
            nothing
        else
            data = split(line[3:10], " ")
            print("data[1] =",data[1],"   data[2] =",data[2])
            push!(X_coor, float(data[1]))
            push!(Y_coor, float(data[2]))
        end

        line_number = line_number + 1
    end
end



body_pts = [[X_coor[i], Y_coor[i], 0.0] for i in 1:NPTS]
t_hats = [[X_tHat[i], Y_tHat[i], 0.0] for i in 1:NPTS]
n_hats = [[X_nHat[i], Y_nHat[i], 0.0] for i in 1:NPTS]
body = VSB.Boundary(body_pts,t_hats,n_hats)

# === Problem parameters ===
println("   CALCULATING PROBLEM PARAMETERS...")
magU_inf = 1.0
U_inf(X) = magU_inf .* [1,0,0]
params = VSB.Parameters(body, U_inf, verbose = true)
alphas = params.alphas

# === Numerical solution for gamma ===
println("   CALCULATING GAMMA NUMERICAL...")
function gamma_numerical_i(i)
    point = body.bodyPTS[i]
    gamma = VSB.CalcVS(body, alphas, point)
    return gamma
end
gammaN = [gamma_numerical_i(i) for i in 1:NPTS]

# === Plotting ===
println("   PLOTTING...")
plot!(theta, gammaN, line = (:blue, :dash), label = "NACA")
savefig(plt1, "VSB_Verification_NACA.pdf")


println("=== END OF Verifications5.jl ===")
