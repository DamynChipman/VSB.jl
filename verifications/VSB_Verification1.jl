# ===== PM2D Verifications: Pressure Distribution =====

println("=== BEGINNING Verifications1.jl ===")

# == Imports ==
using VSB
using SciTools
using LinearAlgebra
using Plots
pyplot()
using LaTeXStrings

cd("/Users/Damyn/Documents/BYU/FLOW Lab/VSB/verifications")

# == Set up geometry and grids ==
N = 100                             # Number of points
NPAN = N - 1                        # Number of panels
del_theta = (2*pi)/(N-1)            # Spacing in theta
R = 1.0                             # Radius
X = zeros(N)                        # X-coordinates
Y = zeros(N)                        # Y-coordinates
theta = zeros(NPAN-1)               # Angle from positive X-axis
CP_analytical = zeros(NPAN-1)       # Analytical CP
CP_PM2D = zeros(NPAN-1)             # PM2D (Semi-Analytical) CP
CP_SVPM = zeros(NPAN-1)             # SVPM CP
for i=1:NPAN
    X[i] = R*cos((i-1)*del_theta)
    Y[i] = R*sin((i-1)*del_theta)
    #println("X[i] = ",X[i],"   Y[i] = ",Y[i])
    if i != N-1
        theta[i] = (i-1)*del_theta
    end
end
X[N] = 1.0; Y[N] = 0.0

# == CP: Analytical ==
for i=1:NPAN-1
    CP_analytical[i] = 1 - 4*sin(theta[i])^2
end

# == CP: VSB w/ Summation of Linear Vortex Lines (Semi-Analytical) ==
oper_cond = [10.0,0.0]
panels = SciTools.Body2D(X,Y)
RHS = zeros(NPAN-1)
for i=1:NPAN-1
    RHS[i] = dot(panels[i].t_hat,oper_cond)
end
alpha_coefs = VSB.CalcVortexSheetCoef(panels,RHS)

# Wrapper function for
function u_gamma(X)
    U = VSB.CalcVortexSheetVelocity(panels,alpha_coefs,X,sigma=1.0,a=10.0)
    return U
end

for i=1:NPAN-1
    CP_PM2D[i] = 1 - (norm(u_gamma([X[i],Y[i],0.0]))/(oper_cond[1]))^2
end

# == CP: PM2D w/ SimpleVPM ==

# Load SimpleVPM
using SimpleVPM
vpm = SimpleVPM


file_path = splitdir(@__FILE__)[1]      # Path to this file


# Output options
save_path = joinpath(file_path, "temps/test00/") # Path where to save paraview files
run_name = "vortex"                     # Name of files
paraview = true                         # Whether to call paraview
prompt = true                           # Whether to prompt the user
verbose = true                          # Whether to print verbose

# ------------- SIMULATION PARAMETERS -------------------------------------
sigma = 1.0                # Smoothing radius
nu = 1e1                  # Kinematic viscosity

magUinf = 10.0
Uinf(t) = magUinf*[1, 0, 0]

# -------------- PARTICLE FIELD --------------------------------------------
# Initiates the field
pfield = ParticleField(nu,
                       zeta_gauserf,       # Basis function
                       g_gauserf;          # Regularizing function of zeta
                       transposed=true,    # Transposed scheme
                       relax=!true,        # Pedrizzeti relaxation scheme
                       rlxf=0.3,           # Relaxation factor
                       integration="euler",# Time integration scheme
                       Uinf=Uinf)


# Add particles on geometry of circle
for i=1:NPAN-1
    loc = [X[i], Y[i], 0.0]
    Gamma = [0.0, 0.0, VSB.CalcVortexSheet(panels,alpha_coefs,loc[1:2])]
    vol = 1.0
    particle = Particle(loc, Gamma, sigma, vol)
    addparticle(pfield, particle)
end

# Calculate pressure coefficient using SVPM velocity function
for i=1:NPAN-1
    CP_SVPM[i] = 1 - (norm(U(pfield, [X[i],Y[i],0.0]))/(oper_cond[1]))^2
end

# == Plotting Results ==
println("=== BEGINNING PLOTTING ===")

plt1 = plot(theta,CP_analytical,line = (:red),
                                label = "Analytical")
plot!(theta,CP_PM2D,line = (:blue, :dash),
                    label = "Semi-Analytical")
plot!(theta,-CP_SVPM,line = (:green, :dash),
                     title = "VSB Verification - Pressure Distribution",
                     xlabel = L"\theta",
                     xlims = (0.0, 2*pi),
                     xticks = ([0.0,pi/2,pi,3*pi/2,2*pi],["0",L"\frac{\pi}{2}",L"\pi",L"\frac{3\pi}{2}",L"2\pi"]),
                     ylabel = L"C_p",
                     ylims = (-4,2),
                     yflip = true,
                     label = "Coupled")
display(plt1)
savefig(plt1,"VSB_Verification_CP.pdf")
