# ===== VSB Verifications: Vortex Sheet Calculations =====

println("=== BEGINNING Verifications5.jl ===")

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
NPTS = 50                              # Number of Points on Boundary
del_theta = (2*pi)/(NPTS - 1)           # Step in theta around body
R = 1.0                                 # Cylinder Radius
X_coor = zeros(NPTS)                    # X-coordinates
Y_coor = zeros(NPTS)                    # Y-coordinates
X_tHat = zeros(NPTS)                    # X Tangent Vector
Y_tHat = zeros(NPTS)                    # X Tangent Vector
X_nHat = zeros(NPTS)                    # Y Normal Vector
Y_nHat = zeros(NPTS)                    # Y Normal Vector
theta = zeros(NPTS)                     # Angle from positive X-axis
for i=1:NPTS
    X_coor[i] = R*cos((i-1)*del_theta)
    Y_coor[i] = R*sin((i-1)*del_theta)

    X_tHat[i] = sin((i-1)*del_theta)
    Y_tHat[i] = -cos((i-1)*del_theta)

    X_nHat[i] = cos((i-1)*del_theta)
    Y_nHat[i] = sin((i-1)*del_theta)

    theta[i] = (i-1)*del_theta
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

# === Analytical solution for gamma ===
println("   CALCULATING GAMMA ANALYTICAL...")
function gamma_analytical_i(i)
    point = body.bodyPTS[i]
    t_hat = body.tHats[i]
    gamma = -dot(U_inf(point), t_hat)
    return gamma
end
gammaA = [gamma_analytical_i(i) for i in 1:NPTS]

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
plt1 = plot(theta, gammaA, line = (:red), label = "Analytical")
plot!(theta, gammaN, line = (:blue, :dash), label = "Numerical")
plot!(xlims = (0.0, 2*pi),
      xticks = ([0.0,pi/2,pi,3*pi/2,2*pi],["0",L"\frac{\pi}{2}",L"\pi",L"\frac{3\pi}{2}",L"2\pi"]),
      xlabel = L"\theta",
      ylims = (-1.05, 1.05),
      ylabel = L"\gamma(\theta)",
      grid = false)
# annotate!( [ (1,0.1,text("Analytical",14,:red)) , (1,-0.1,text("Numerical",14,:blue))])
savefig(plt1, "VSB_Verification_VS2.pdf")


println("=== END OF Verifications5.jl ===")
