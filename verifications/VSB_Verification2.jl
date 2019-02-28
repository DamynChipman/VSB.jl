# ===== VSB Verifications: Slip Velocity =====

println("=== BEGINNING Verifications2.jl ===")

# === Imports ===
println("   IMPORTING PACKAGES...")
using VSB
using SimpleVPM
using LinearAlgebra
using Printf

cd("/Users/Damyn/Documents/BYU/FLOW Lab/VSB/verifications")

# === Geometry and Grids ===
println("   GENERATING GEOMETRY...")
NPTS = 500                              # Number of Points on Boundary
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

# === RBF Coefficients and Velocity Functions ===
U_inf = [1.0, 0.0, 0.0]
U_inf_func(X) = U_inf
print("   CALCULATING ALPHAS...")
alphas = VSB.CalcVSCoefs(body, U_slip=U_inf_func)
print("   DONE!")
println()

U_gamma(X) = VSB.CalcVSVelocityCirlce(body, alphas, X)

println("   CALCULATING VELOCITIES...")
println(" -      X     -  |  - U_SLIP TANGENT -  |  - U_SLIP NORMAL -  |  -      U_GAMMA      -  |  - U_GAMMA TANGENT -  |  - U_GAMMA NORMAL -  |")
println("=========================================================================================================================================")
# === Calculate U_slip along surface ===
for i in 1:NPTS
    X_eval = body.bodyPTS[i]
    n_hat = body.nHats[i]
    t_hat = body.tHats[i]

    u_gamma = U_gamma(X_eval)

    u_slip_normal = dot(U_inf_func(X_eval), n_hat)
    u_gamma_normal = dot(U_gamma(X_eval), n_hat)

    u_slip_tangent = dot(U_inf_func(X_eval), t_hat)
    u_gamma_tangent = dot(U_gamma(X_eval), t_hat)

    if i % 10 == 0
        println("  [",@sprintf("%5.2f",X_eval[1]),",",@sprintf("%5.2f",X_eval[2]),"]",
                "  |  ",@sprintf("%18.3f",u_slip_tangent),
                "  |  ",@sprintf("%17.3f",u_slip_normal),
                "  |  [",@sprintf("%5.2f",u_gamma[1]),", ",@sprintf("%5.2f",u_gamma[2]),", ",@sprintf("%5.2f",u_gamma[3]),"]",
                "  |  ",@sprintf("%19.3f",u_gamma_tangent),
                "  |  ",@sprintf("%18.3f",u_gamma_normal),
                "  |")
    end
end
println("=========================================================================================================================================")

println("=== END OF Verifications2.jl ===")
