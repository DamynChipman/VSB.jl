# ===== PM2D Verifications: Pressure Distribution =====
println("HERE")
# == Imports ==
using VSB
using SciTools
println("HERE")
# == Set up geometry and grids ==
N = 200                             # Number of points
NPAN = N - 1                        # Number of panels
del_theta = (2*pi)/(N-1)            # Spacing in theta
R = 1.0                             # Radius
X = zeros(N)                        # X-coordinates
Y = zeros(N)                        # Y-coordinates
theta = zeros(NPAN-1)               # Angle from positive X-axis
CP_analytical = zeros(NPAN-1)       # Analytical CP
CP_PM2D = zeros(NPAN-1)             # PM2D (Semi-Analytical) CP
CP_SVPM = zeros(NPAN-1)             # SVPM CP
for i=1:NPAN-1
    X[i] = R*cos((i-1)*del_theta)
    Y[i] = R*sin((i-1)*del_theta)
    if i != N-1
        theta[i] = (i-1)*del_theta
    end
end

# == CP: Analytical ==
for i=1:NPAN-1
    CP_analytical[i] = 1 - 4*sin(theta[i])^2
end

# == CP: VSB w/ Summation of Linear Vortex Lines (Semi-Analytical) ==
oper_cond = [10.0,0.0]
panels = SciTools.Body2D(X,Y)
RHS = ones(NPAN) * 10.0
alpha_coefs = VSB.CalcVortexSheetCoef(panels,RHS)

# Wrapper function for
function u_gamma(X)
    U = VSB.CalcVortexSheetVelocity(panels,alpha_coefs,X,sigma=0.2,a=18.5*pi)
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
sigma = 1.0               # Smoothing radius
nu = 1e1                  # Kinematic viscosity

magUinf = 0
Uinf(t) = magUinf*[1, 0, 0]
nsteps = 150              # Number of time steps
dt = 0.001                # Time step

nsteps_relax = 1          # Steps in between relaxation
sgm0 = sigma              # Default core size
beta_cs = 1.1            # Maximum core size growth

# -------------- PARTICLE FIELD --------------------------------------------
# Initiates the field
pfield = ParticleField(nu,
                        zeta_gauserf,       # Basis function
                        g_gauserf;          # Regularizing function of zeta
                        transposed=true,    # Transposed scheme
                        relax=!true,        # Pedrizzeti relaxation scheme
                        rlxf=0.3,           # Relaxation factor
                        integration="euler",# Time integration scheme
                        Uinf=Uinf
                      )


for i=1:NPAN-1
    loc = [X[i], Y[i], 0.0]
    Gamma = alpha_coefs[i]
    vol = 1.0
    particle = Particle(loc, Gamma, sigma, vol)
    addparticle(pfield, particle)
end

for i=1:NPAN-1
    CP_SVPM[i] = 1 - (norm(U(pfield, [X[i],Y[i],0.0]))/(oper_cond[1]))^2
end
