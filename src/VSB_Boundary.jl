"""
`Boundary`


"""
mutable struct Boundary

    bodyPTS::NUMB_MATRIX
    tHats::NUMB_MATRIX
    nHats::NUMB_MATRIX
    NPTS_BODY::Int

    function Boundary(body_pts::NUMB_MATRIX,
                      t_hats::NUMB_MATRIX,
                      n_hats::NUMB_MATRIX)
        new(body_pts,t_hats,n_hats,length(body_pts))
    end
end

"""

"""
function CalcRhoCoefs(self::Boundary)

    # === Extract geometry from boundary ===
    X_body = [[point[1], point[2], 0.0] for point in self.bodyPTS]
    n_hats = [[n_hat[1], n_hat[2], 0.0] for n_hat in self.nHats]
    t_hats = [[t_hat[1], t_hat[2], 0.0] for t_hat in self.tHats]
    NPTS = self.NPTS_BODY

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

    # === Coefficient functions ===
    function Theta_ki(k,i)
        res = 0
        for j in 1:NPTS
            if j != k
                res = res + (1/(2*pi)) * (dot(R_kj(k,j), n_hats[j]) * RBF_gauss(r_ji(j,i)) * del_S(j))/(r_kj(k,j)^2)
            end
        end
        return res
    end

    function phi_ki(k,i)
        return RBF_gauss(r_ki(k,i))
    end

    # === Build coefficient matrix ===
    A = zeros(NPTS,NPTS)
    for k in 1:NPTS
        for i in 1:NPTS
            A[k,i] = Theta_ki(k,i) - phi_ki(k,i)
        end
    end

    # === Build RHS vector ===
    b = zeros(NPTS)

    # === Solve linear system ===
    etas = A\b
    return etas

end

"""

"""
function CalcRho(self::Boundary, etas::Array{T}, X_eval::Array{T}) where {T<:Real}

    # === Unpack geometry ===
    XP = [[point[1], point[2], 0.0] for point in boundary.bodyPTS]
    NPTS = boundary.NPTS_BODY

    # === Summation over all points for gamma ===
    rho = sum( [etas[i] * RBF_gauss(norm(X_eval - XP[i])) for i in 1:NPTS ] )
    return rho

end
