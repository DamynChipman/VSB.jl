# ===== VSB Verifications: Velocity Profiles =====

println("=== BEGINNING Verifications4.jl ===")

# === Imports ===
println("   IMPORTING PACKAGES...")
using VSB
using Plots
pyplot()

cd("/Users/Damyn/Documents/BYU/FLOW Lab/VSB/verifications")

# === Geomerty of problem ===
println("   GENERATING GEOMETRY...")
NPTS = 250                              # Number of Points on Boundary
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

# === Vortex sheet calculations ===
println("   CALCULATING PROBLEM PARAMETERS...")
magU_inf = 1.0
U_inf(X) = magU_inf .* [1.0, 0.0, 0.0]

params = VSB.Parameters(body, U_inf)
alphas = params.alphas
U(X) = VSB.CalcVSVelocityCirlce(body, alphas, X)

# === Generate velocity field ===
println("   CALCULATING VELOCITY FIELD...")



println("=== END OF Verifications1.jl ===")
