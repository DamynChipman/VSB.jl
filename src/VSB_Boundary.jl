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

    println(" NACA digits: m = ",m,"  p = ",p,"  tau = ",tau)

    # Define camber geometery
    if p != 0
        del_x1 = (p*c - 0)/(N/4)
        del_x2 = (c - p*c)/(N/4)

        x1 = 0:del_x1:(p*c)
        Y1 = ((m)./(p^2)).*(2*p*(x1./c) - (x1./c).^2)
        dY_dx1 = ((2*m)/(p^2)).*(p - x1./c)

        x2 = (p*c):del_x2:(c)
        Y2 = ((m)./((1 - p)^2)).*((1 - 2*p) + (2*p).*(x2./c) - (x2./c).^2)
        dY_dx2 = ((m)/((1 - p)^2)).*((1 - 2*p) + 2*p.*(x2./c) - (x2./c).^2)

        camber = [x1 Y1; x2 Y2]

        x = zeros(length(x1)+length(x2)-1)
        Y = zeros(length(Y1)+length(Y2)-1)
        dY_dx = zeros(length(dY_dx1)+length(dY_dx2)-1)
        for i=1:length(x1)
            x[i] = x1[i]
            Y[i] = Y1[i]
            dY_dx[i] = dY_dx1[i]
        end
        for i=2:length(x2)
            x[i+length(x1)-1] = x2[i]
            Y[i+length(x1)-1] = Y2[i]
            dY_dx[i+length(x1)-1] = dY_dx2[i]
        end
    else
        del_x = c/(N/2)
        x = 0:del_x:c
        Y = zeros(length(x),1)
        camber = [x Y]
        dY_dx = zeros(length(x),1)
    end

    # Define airfoil thickness
    T(x) = 5*tau*(0.2969*x^.5 - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)

    # Determine x and y coordinates of airfoil surface
    x_upper, x_lower, y_upper, y_lower = 0. * x, 0. * x, 0. * x, 0. * x
    for i=1:length(x)
        theta = atan(dY_dx[i])

        x_upper[i] = x[i] - T(x[i])*sin(theta)
        x_lower[i] = x[i] + T(x[i])*sin(theta)
        y_upper[i] = Y[i] + T(x[i])*cos(theta)
        y_lower[i] = Y[i] - T(x[i])*cos(theta)
    end

    upper = zeros(length(x_upper),2)
    lower = zeros(length(x_lower),2)
    for i=1:length(upper[:,1])
        upper[i,:] = [x_upper[i],y_upper[i]]
        lower[i,:] = [x_lower[i],y_lower[i]]
    end

    if m == 0 && p == 0
        upper[1,:] = [0.0 0.0]
        lower[1,:] = [0.0 0.0]
    end

    # Redefine the orientation of the airfoil
    airfoil = zeros(2*length(x_lower)-2,2)
    for n=1:length(upper[:,1])
        airfoil[n,1] = upper[(length(upper[:,1])+1)-n,1]
        airfoil[n,2] = upper[(length(upper[:,1])+1)-n,2]
    end
    for n=2:(length(upper[:,1])-1)
        airfoil[n+length(upper[:,1])-1,1] = lower[n,1]
        airfoil[n+length(upper[:,1])-1,2] = lower[n,2]
    end

    return upper,lower,airfoil,camber
end
