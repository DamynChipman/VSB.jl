# ===== VSB Verifications: RunTime_Function Building =====

println(" === BEGINNING Verifications3.jl === ")

# === Imports ===
println("   IMPORTING PACKAGES...")
using VSB
using SimpleVPM
using LinearAlgebra
using Printf

cd("/Users/Damyn/Documents/BYU/FLOW Lab/VSB/verifications")

# === Boundary Geometry ===
println("   GENERATING GEOMETRY...")
NPTS = 10                               # Number of Points on Boundary
del_theta = (2*pi)/(NPTS    )           # Step in theta around body
R = 1.0                                 # Cylinder Radius
k_hat = [0.0, 0.0, 1.0]                 # Z direction unit vector
X_coor = zeros(NPTS)                    # X-coordinates
Y_coor = zeros(NPTS)                    # Y-coordinates
X_tHat = zeros(NPTS)                    # X Tangent Vector
Y_tHat = zeros(NPTS)                    # X Tangent Vector
X_nHat = zeros(NPTS)                    # Y Normal Vector
Y_nHat = zeros(NPTS)                    # Y Normal Vector
overlap = 2.0                           # Particle overlap
sigma = overlap * del_theta             # Guassian Spreading
phase = 0.0 * del_theta
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

# === Parameters ===
# Run 5
nu = 1e0                                    # Kinematic Viscosity
magU_inf = 1.0                              # Magnitude of Free-Stream Velocity
del_t = 0.5                                 # Time step
NSTEPS = 100                                # Number of time steps
Re = magU_inf * R / nu                      # Reynolds Number

U_inf(t) = magU_inf * [1.0, 0.0, 0.0]       # Free-Stream velocity

# === Set up SimpleVPM ===
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

# === Initial vorticity field ===
params = VSB.Parameters(body, U_slip, verbose = false, sigma=sigma)
alphas = params.alphas
betas = VSB.CalcVSDiffusionCoefs(body, alphas, del_t, pfield.nu, sigma=sigma)
betas[1] = 0.0
betas[NPTS] = 0.0

offset = 0.01             # Normal direction offset for new particle locations
for i in 1:NPTS
    loc = body_pts[i] + offset .* n_hats[i]
    Gamma = (1/sigma^3) * betas[i] .* k_hat
    vol = 1.0
    particle = SimpleVPM.Particle(loc, Gamma, sigma, vol)
    SimpleVPM.addparticle(pfield, particle)
end

# === RunTimeFunction ===
function RunTimeFunction(pfield, t, del_t)

    if (t != 0)

        for p in pfield.particles
            p.Gamma[1] = 0.0
            p.Gamma[2] = 0.0
        end

        # === Calculate parameters ===
        params = VSB.Parameters(body, U_slip, verbose = false)
        alphas = params.alphas
        betas = VSB.CalcVSDiffusionCoefs(body, alphas, del_t, pfield.nu)

        # === Add particles on surface ===
        offset = 0             # Normal direction offset for new particle locations
        for i in 1:NPTS
            loc = body_pts[i] + offset .* n_hats[i]
            Gamma = (1/sigma^3) * betas[i] .* k_hat
            vol = 1.0
            particle = SimpleVPM.Particle(loc, Gamma, sigma, vol)
            SimpleVPM.addparticle(pfield, particle)
        end

    end

    # === Returns ===
    return false

end

# === Simulation Parameters ===
println("   RUNNING SIMULATION...")

println("***********************************")
println(" ===== SIMULATION PARAMETERS ===== ")
println(" NU:      ",nu)
println(" U_INF:   ",magU_inf)
println(" NSTEPS:  ",NSTEPS)
println(" NPTS:    ",NPTS)
println(" DEL_T:   ",del_t)
println(" RE:      ",Re)
println("***********************************")

NSTEPS_relax = -1
sigma0 = sigma
file_path = splitdir(@__FILE__)[1]
save_path = joinpath(file_path, "sims/verification3_07/")
run_name = "verifications3_07"
paraview = true

SimpleVPM.run_vpm!(pfield, del_t, NSTEPS,
                   runtime_function = RunTimeFunction,
                   nsteps_relax = NSTEPS_relax,
                   rbf_itmax = 1000,
                   sgm0 = sigma0,
                   save_path = save_path,
                   run_name = run_name)

if paraview
    println("   CALLING PARAVIEW...")
    strn = run_name * "...xmf;"
    run(`paraview --data="$save_path$strn"`)
end

println(" === END OF Verifications3.jl === ")
