"""
`CalcVSCoefs(boundary, etas; U_slip, sigma)`


"""
function CalcVSCoefs(boundary::Boundary,
                     etas::Array{T};
                     U_slip::Union{Nothing, Function}=nothing,
                     sigma=0.2) where {T<:Real}

    # === Extract geometry from boundary ===
    X_body = [[point[1], point[2], 0.0] for point in boundary.bodyPTS]
    n_hats = [[n_hat[1], n_hat[2], 0.0] for n_hat in boundary.nHats]
    t_hats = [[t_hat[1], t_hat[2], 0.0] for t_hat in boundary.tHats]
    NPTS = boundary.NPTS_BODY

    # === Function definitions ===

    # i -> RBF index
    # j -> Integration index
    # k -> Evaluation index

    X_i(i) = X_body[i]                                  # Set_i of points on S
    X_j(j) = X_body[j]                                  # Set_j of points on S
    X_jp1(j) = (j != NPTS) ? X_body[j+1] : X_body[1]    # Point j + 1
    X_k(k) = X_body[k]                                  # Set_k of points on S

    R_kj(k,j) = X_k(k) - X_j(j)                         # X_k - X_j
    R_ji(j,i) = X_j(j) - X_i(i)                         # X_j - X_i
    R_ki(k,i) = X_k(k) - X_i(i)                         # X_k - X_i
    r_kj(k,j) = norm(R_kj(k,j))                         # | X_k - X_j |
    r_ji(j,i) = norm(R_ji(j,i))                         # | X_j - X_i |
    r_ki(k,i) = norm(R_ki(k,i))                         # | X_k - X_i |

    del_S(j) = norm(X_j(j) - X_jp1(j))                  # ΔS_j
    L = sum( [del_S(j) for j in 1:NPTS] )               # Length of surface

    #etas = CalcRhoCoefs(boundary)                       # Coefficients for rho RBF
    rho_1(k) = CalcRho(boundary, etas, X_k(k))          # ρ_1 of Eigendecomposition

    # === Coefficent functions ===
    function Theta_ki(k,i)
        res = 0
        for j in 1:NPTS
            if j != k
                res = res + (1/(2*pi)) * (dot(R_kj(k,j), n_hats[j]) * RBF_gauss(r_ji(j,i), sigma=sigma) * del_S(j))/(r_kj(k,j)^2)
            end
        end
        return res
    end
    #Theta_ki(k,i) = sum([((k != j) ? (1/pi) * (dot(R_kj(k,j), n_hats[j]) * RBF_gauss(r_ji(j,i)) * del_S(j))/(r_kj(k,j)^2) : 0) for j in 1:NPTS] )
    Lambda_ki(k,i) = sum( [RBF_gauss(r_ji(j,i), sigma=sigma) * (rho_1(k))/(L) for j in 1:NPTS] )
    phi_ki(k,i) = RBF_gauss(r_ki(k,i), sigma=sigma)

    # === Calculate coefficient matrix ===
    A = zeros(NPTS,NPTS)
    for k in 1:NPTS
        for i in 1:NPTS
            A[k,i] = phi_ki(k,i) - Theta_ki(k,i) + Lambda_ki(k,i)
        end
    end
    #A = [ [phi_ki(k,i) - Theta_ki(k,i) + Lambda_ki(k,i) for k in 1:NPTS] for i in 1:NPTS]
    #println("A: ",A)

    # === Calculate RHS vector ===
    b = zeros(NPTS)
    for i in 1:NPTS
        b[i] = dot(U_slip(X_i(i)), t_hats[i])
    end
    #b = [dot(U_slip(X_i(i)), t_hats[i]) for i in 1:NPTS]
    #println("b: ",b)

    # === Calculate alpha coefficients ===
    alphas = A\b
    return alphas

end

"""
`CalcVortexSheet(panels,alpha,X,sigma=0.2,a=1.0)`

Given a geometry of panels, the strength of the vortex particle on the panels,
and an X location, computes the contribution from all the vortex particles on the
surface at the given X location.

# ARGUMENTS
* `panels::Array{Panel2D}` : Array of Panel2D objects representing a body
* `alpha::Array{Float64`   : Radial basis function coefficients previously solved for
* `X::Array{Float64}`      : Point to compute contribution at
* `sigma::Float64=0.2`     : Guassian width
* `a::Float64=1.0`         : Normalization constant for Guassian

# OUTPUTS
* `gamma::Float64`         : Contribution from all vortex particles at X
"""
function CalcVS(boundary::Boundary,
                alpha::Array{T},
                X_eval::Array{T};
                sigma=0.2) where {T<:Real}

    # === Unpack geometry ===
    XP = [[point[1], point[2], 0.0] for point in boundary.bodyPTS]
    NPTS = boundary.NPTS_BODY

    # === Summation over all points for gamma ===
    gamma = -sum( [alpha[i] * RBF_gauss(norm(X_eval - XP[i]), sigma=sigma) for i in 1:NPTS ] )
    return gamma

end

"""

"""
function CalcVSVelocityCirlce(boundary::Boundary,
                              alpha::Array{T},
                              X_eval::Array{T};
                              radius=1.0) where {T<:Real}

    # === Unpack geometry ===
    X_body = [[point[1], point[2], 0.0] for point in boundary.bodyPTS]
    NPTS = boundary.NPTS_BODY

    # === Function definitions ===
    X_j(j) = X_body[j]                                       # Set_j of points on S
    X_jp1(j) = (j != NPTS) ? X_body[j+1] : X_body[1]         # Point j + 1
    R_xj(j) = X_eval - X_j(j)                                # X_eval - X_j
    r_xj(j) = norm(R_xj(j))                                  # | X_eval - X_j |

    gamma(j) = CalcVS(boundary, alpha, X_j(j)) .* [0, 0, 1]  # Vortex sheet strength
    del_S(j) = norm(X_j(j) - X_jp1(j))                       # ΔS_j

    num(j) = cross(X_eval - X_j(j), gamma(j)) * del_S(j)     # Numerator
    den(j) = 2*pi * r_xj(j)^2                                # Denominator

    # === Calculate velocity ===
    U = sum( [((r_xj(j) < 1e-12) ? zeros(3) : num(j) ./ den(j)) for j in 1:NPTS] )
    return U

end
