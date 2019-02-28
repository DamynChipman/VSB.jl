# ===== VSB Verifications: Pressure Distribution =====

println("=== BEGINNING Verifications1.jl ===")

# == Imports ==
println("   IMPORTING PACKAGES...")
using VSB
using SimpleVPM
using LinearAlgebra
using Plots
pyplot()
using LaTeXStrings
using Printf

cd("/Users/Damyn/Documents/BYU/FLOW Lab/VSB/verifications")

# === Boundary Geometry and Plotting Grids ===
println("   GENERATING GEOMETRY...")
NPTS = 150                              # Number of Points on Boundary
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

# === CP: Analytical ===
println("   CALCULATING CP_ANALYTICAL...")
CP_analytical = [1 - 4*sin(theta[i])^2 for i in 1:NPTS]

# === CP: VSB ===
println("   CALCULATING CP_VSB...")
println("   CALCULATING PROBLEM PARAMETERS...")

U_inf = [1.0, 0.0, 0.0]
U_inf_func(X) = U_inf

params = VSB.Parameters(body, U_inf_func, verbose = true)
alphas = params.alphas

U_gamma(X) = VSB.CalcVSVelocityCirlce(body, alphas, X)

CP_VSB = zeros(NPTS)
for i in 1:NPTS
    if i % 10 == 0
        println("     - ON POINT ",i," OF ",NPTS)
    end
    scale_factor = 0.000
    X_eval = body.bodyPTS[i] + scale_factor .* body.nHats[i]
    CP_VSB[i] = 1 - (norm( U_gamma(X_eval)) / norm( U_inf_func(X_eval) ))^2
end

# === CP: VPM ===
println("   CALCULATING CP_VPM...")
CP_VPM = zeros(NPTS)

nu = 1e1
Uinf(t) = U_inf
sigma = 1.0

print("     - CREATING PARTICLE FIELD...")
pfield = SimpleVPM.ParticleField(nu,                 # Kinematic viscosity
                                 zeta_gauserf,       # Basis function
                                 g_gauserf;          # Regularizing function of zeta
                                 transposed=false,    # Transposed scheme
                                 relax=!true,        # Pedrizzeti relaxation scheme
                                 rlxf=0.3,           # Relaxation factor
                                 integration="euler",# Time integration scheme
                                 Uinf=Uinf)          # Free stream velocity
print("   DONE!")
println()

print("     - ADDING PARTICLES...")
for i=1:NPTS
    loc = body.bodyPTS[i]
    Gamma = [0.0, 0.0, VSB.CalcVS(body,alphas,loc)]
    vol = 1.0
    particle = Particle(loc, Gamma, sigma, vol)
    addparticle(pfield, particle)
end
print("   DONE!")
println()

for i=1:NPTS
    if i % 10 == 0
        println("     - ON POINT ",i," OF ",NPTS)
    end
    CP_VPM[i] = -(1 - ( norm( U(pfield, body.bodyPTS[i]) ) / norm( U_inf ))^2)
end

println()
println("   RESULTS...")
println("CP_ANALYTICAL    CP_VSB    CP_VPM")
for i in 1:NPTS
    if i % 10 == 0
        println(@sprintf("%13.3f",CP_analytical[i]),"  ",
                @sprintf("%8.3f",CP_VSB[i]),"  ",
                @sprintf("%8.3f",CP_VPM[i]))
    end
end

# === Plotting ===
println("   PLOTTING...")
plt1 = plot(theta, CP_analytical, line = (:red), label = "Analytical")
plot!(theta, CP_VSB, line = (:blue, :dash), label = "VSB")
plot!(title = "VSB Verification - Pressure Distribution",
      xlabel = L"\theta",
      xlims = (0.0, 2*pi),
      xticks = ([0.0,pi/2,pi,3*pi/2,2*pi],["0",L"\frac{\pi}{2}",L"\pi",L"\frac{3\pi}{2}",L"2\pi"]),
      ylabel = L"C_P",
      ylims = (-4,2),
      yflip = true)

plt2 = plot(theta, alphas, line=(:red), title = "Alpha Coefficients")
savefig(plt1, "VSB_Verification_CP.pdf")
savefig(plt2, "VSB_Verification_Alphas.pdf")

println("=== END OF Verifications1.jl ===")
