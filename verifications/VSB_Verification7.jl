# ===== VSB Verifications: Vorticity Coupling of VSB and VPM =====

println(" === BEGINNING Verifications7.jl === ")

# === Imports ===
println("   IMPORTING PACKAGES...")
using VSB
using SimpleVPM
using LinearAlgebra
using Printf

cd("/Users/Damyn/Documents/BYU/FLOW Lab/VSB/verifications")

# === Boundary Geometry ===
println("   GENERATING GEOMETRY...")
NPTS = 20                              # Number of Points on Boundary
del_theta = (2*pi)/(NPTS - 1)           # Step in theta around body
R = 1.0                                 # Cylinder Radius
k_hat = [0.0, 0.0, 1.0]                 # Z direction unit vector
X_coor = zeros(NPTS)                    # X-coordinates
Y_coor = zeros(NPTS)                    # Y-coordinates
X_tHat = zeros(NPTS)                    # X Tangent Vector
Y_tHat = zeros(NPTS)                    # X Tangent Vector
X_nHat = zeros(NPTS)                    # Y Normal Vector
Y_nHat = zeros(NPTS)                    # Y Normal Vector
phase = 0.0 * del_theta                 # Phase shift
for i=1:NPTS
    X_coor[i] = R*cos((i-1)*del_theta + phase)
    Y_coor[i] = R*sin((i-1)*del_theta + phase)

    X_tHat[i] = sin((i-1)*del_theta + phase)
    Y_tHat[i] = -cos((i-1)*del_theta + phase)

    X_nHat[i] = cos((i-1)*del_theta + phase)
    Y_nHat[i] = sin((i-1)*del_theta + phase)
end
body_pts = [[X_coor[i], Y_coor[i], 0.0] for i in 1:NPTS]
t_hats = [[X_tHat[i], Y_tHat[i], 0.0] for i in 1:NPTS]
n_hats = [[X_nHat[i], Y_nHat[i], 0.0] for i in 1:NPTS]
body = VSB.Boundary(body_pts,t_hats,n_hats)

# === Set up SimpleVPM ===
println("   SETTING UP SVPM...")
nu = 1e0                                                      # Kinematic Viscosity
dt = 1e0                                                      # Time step (place holder here)
overlap = 2.0                                                 # Particle overlap
sigma = overlap * del_theta                                   # Guassian Spreading
magU_inf = 1.0                                                # Magnitude of Free-Stream Velocity
U_inf(t) = magU_inf * [1.0, 0.0, 0.0]                         # Free-Stream velocity

println("   CREATING PARTICLE FIELD...")
pfield = SimpleVPM.ParticleField(nu,
                                 zeta_gauserf,                # Basis function
                                 g_gauserf;                   # Regularizing function of zeta
                                 transposed=true,             # Transposed scheme
                                 relax=!true,                 # Pedrizzeti relaxation scheme
                                 rlxf=0.3,                    # Relaxation factor
                                 integration="euler",         # Time integration scheme
                                 Uinf=U_inf)                  # Free-Stream velocity function

# === Velocity function ===
function U_slip(X)
    if length(pfield.particles) == 0
        return U_inf(0)
    else
        return SimpleVPM.U(pfield, X) + U_inf(0)
    end
end

# === Calculate parameters ===
println("   CALCULATING PARAMS...")
params = VSB.Parameters(body, U_slip, verbose = false, sigma=sigma)
alphas = params.alphas
betas = VSB.CalcVSDiffusionCoefs(body, alphas, dt, pfield.nu, sigma=sigma)
betas[1] = 0.0
betas[NPTS] = 0.0

println("betas: ",betas)

# === Add particles on surface ===
println("   ADDING PARTICLES TO FIELD...")
offset = 0             # Normal direction offset for new particle locations
for i in 1:NPTS
    X = body_pts[i] + offset .* n_hats[i]
    Gamma = betas[i] .* k_hat
    #Gamma = 1.0 .* k_hat
    vol = 1.0
    particle = SimpleVPM.Particle(X, Gamma, sigma, vol)
    SimpleVPM.addparticle(pfield, particle)
end

SimpleVPM.save(pfield,"pfield_testing")

# === Vorticity calculations ===
println("   BEGIN VORTICITY CALCULATIONS:")
println()
println()
# Method 1: VSB from diffused vortex sheet
omega_VSB(X) = (1/sigma^3)*VSB.CalcDiffusion(body, betas, X, sigma=sigma) .* k_hat

# Method 2: VPM RBF
omega_RBF(X) = SimpleVPM.omega_approx(pfield, X)

# Method 3: Definition of vorticity via VPM
function omega_VPM(i)
    junk, vorticity = SimpleVPM.vortexstretching(pfield, pfield)
    return vorticity[i]
end

# === Results ===
for i in 1:NPTS
    X = body_pts[i]
    XP = pfield.particles[i].X
    VSB = omega_VSB(X)
    RBF = omega_RBF(X)
    VPM = omega_VPM(i)

    println("=== i: ",@sprintf("%3.0i",i)," ======================================")
    println("            X = [",@sprintf("%3.2f",X[1]),", ",@sprintf("%3.2f",X[2]),"]")
    println("           XP = [",@sprintf("%3.2f",XP[1]),", ",@sprintf("%3.2f",XP[2]),"]")

    # println("     gamma[i]")
    println("      beta[i] = ",betas[i] .* k_hat)
    println("     Gamma[i] = ",pfield.particles[i].Gamma)

    println("    omega_VSB = [",@sprintf("%6.4f",VSB[1]),", ",@sprintf("%6.4f",VSB[2]),", ",@sprintf("%6.4f",VSB[3]),"]")
    println("    omega_RBF = [",@sprintf("%6.4f",RBF[1]),", ",@sprintf("%6.4f",RBF[2]),", ",@sprintf("%6.4f",RBF[3]),"]")
    println("    omega_VPM = [",@sprintf("%6.4f",VPM[1]),", ",@sprintf("%6.4f",VPM[2]),", ",@sprintf("%6.4f",VPM[3]),"]")

    println("  |omega_VSB| = ",@sprintf("%6.2f",norm(omega_VSB(X))))
    println("  |omega_RBF| = ",@sprintf("%6.2f",norm(omega_RBF(X))))
    println("  |omega_VPM| = ",@sprintf("%6.2f",norm(omega_VPM(i))))
    println()
end
