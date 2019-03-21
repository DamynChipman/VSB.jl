"""
    Boundary(body_pts,t_hats,n_hats)

Struct containing boundary information. The user must supply the location of the
body points (body_pts), and the tangent and normal vectors of the body at each
body point.

# ARGUMENTS
* `body_pts`   : Array of the form [[x1,y1,z1], [x2,y2,z2], ...] containing the
                 locations of the boundary nodes. Assumes a CCW orientation.
* `t_hats`     : Array of the same form containing the tangent unit vectors at
                 each ponit of the boundary. Same size as body_pts.
* `n_hats`     : Array of the same form containing the normal unit vectors at
                 each ponit of the boundary. Same size as body_pts.

# PROPERTIES
* `bodyPTS`    : Array of the form [[x1,y1,z1], [x2,y2,z2], ...] containing the
                 locations of the boundary nodes.
* `tHats`      : Array of the same form containing the tangent unit vectors at
                 each ponit of the boundary. Same size as body_pts.
* `nHats`      : Array of the same form containing the normal unit vectors at
                 each ponit of the boundary. Same size as body_pts.
* `NPTS_BODY`  : Number of points in the body.
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
function CalcRhoCoefs(self::Boundary; sigma = 0.2)

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
                res = res + (1/(2*pi)) * (dot(R_kj(k,j), n_hats[j]) * RBF_gauss(r_ji(j,i), sigma=sigma) * del_S(j))/(r_kj(k,j)^2)
            end
        end
        return res
    end

    function phi_ki(k,i)
        return RBF_gauss(r_ki(k,i), sigma=sigma)
    end

    # === Build coefficient matrix ===
    A = zeros(NPTS,NPTS)
    for k in 1:NPTS
        for i in 1:NPTS
            A[k,i] = Theta_ki(k,i) - phi_ki(k,i)
        end
    end

    #println(A)

    # === Build RHS vector ===
    b = zeros(NPTS)

    # === Solve linear system ===
    etas = A\b
    return etas

end

"""

"""
function CalcRho(self::Boundary, etas::Array{T}, X_eval::Array{T}; sigma = 0.2) where {T<:Real}

    # === Unpack geometry ===
    XP = [[point[1], point[2], 0.0] for point in self.bodyPTS]
    NPTS = self.NPTS_BODY

    # === Summation over all points for gamma ===
    rho = sum( [etas[i] * RBF_gauss(norm(X_eval - XP[i]), sigma=sigma) for i in 1:NPTS ] )
    return rho

end

"""
    `NACA4(numb, N, c)`

NACA Four-Digit Series Airfoil. Returns two lists of ordered pairs representing
a NACA Four-Digit Airfoil.

# ARGUMENTS
* `numb::String`     : Four digit series. Symmetric airfoil given by "00xx"
* `N::Int64`         : Number of points
* `c::Float64=1.0`   : Length of airfoil. Defaults to 1.0 for appropiate scaling
"""
function NACA4(numb::String,
               N::Int64,
               c::Float64=1.0)

    # Unpackage the digits
    m = parse(Int,numb[1])*.01
    p = parse(Int,numb[2])*.1
    tau = (parse(Int,numb[3])*10 + parse(Int,numb[4]))*.01

    Yt(x) = 5*tau*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^5)

    function Yc(x)
        if 0 <= x && x > p*c
            return (m/p^2) * (2*p*(x/c) - (x/c)^2)
        else
            return (m/(1 - p)^2) * ((1 - 2*p) + 2*p*(x/c) - (x/c)^2)
        end
    end

    function dYc_dX(x)
        if 0 <= x && x > p*c
            return (2*m/p^2) * (p - x/c)
        else
            return (2*m/(1 - p)^2) * (p - x/c)
        end
    end

    theta(x) = atan(dYc_dX(x))

    del_X = c/(N - 1)
    X = 0:del_X:c
    X_upper = [x - Yt(x)*sin(theta(x)) for x in X]
    X_lower = [x + Yt(x)*sin(theta(x)) for x in X]
    Y_upper = [Yc(x) + Yt(x)*cos(theta(x)) for x in X]
    Y_lower = [Yc(x) - Yt(x)*cos(theta(x)) for x in X]

    upper = [[X_upper[i], Y_upper[i]] for i in 1:N]
    lower = [[X_lower[i], Y_lower[i]] for i in 1:N]

    body_pts = [[0.0,0.0,0.0] for i in 1:2*N]
    for i in 0:N-1
        body_pts[i+1] = [X_upper[N-i], Y_upper[N-i], 0.0]
        body_pts[N+i+1] = [X_lower[i+1], Y_lower[i+1], 0.0]
    end

    t_hats = [[0.0,0.0,0.0] for i in 1:2*N]
    n_hats = t_hats
    for i in 1:2*N

        if i != 2*N
            R = body_pts[i] - body_pts[i+1]
        else
            R = body_pts[N] - body_pts[1]
        end

        r = norm(R)
        t_hats[i] = R/r
        n_hats[i] = cross(t_hats[i], [0,0,1])
    end

    NACA = Boundary(body_pts, t_hats, n_hats)

    return NACA

end
