"""
`CalcVortexSheetCoef(panels,RHS)`

Given a geometry of LineSegments, computes the strength of a vortex sheet represented
by particles collcated on the starting and ending points of the panel. Builds
the coefficient matrix and performs the numerical integration with Julia's QuadGK
quadrature package. Returns a vector of radial basis function coefficients for
the strength of the vortex sheet required to induce a no-slip velocity on the body.

# ARGUMENTS
* `panels::Array{LineSegment}` : Array of Panel2D objects representing a body
* `RHS::Array{Float64}`        : Given RHS vector for slip velocity

# OUTPUTS
* `alpha::Array{Float64}`  : Solution vector of coefficients for RBF vortex sheet strength
"""
# function CalcVortexSheetCoef(panels::Array{SciTools.LineSegment},
#                              RHS::Array{Float64})
function CalcVSCoef(boundary::Boundary,
                    RHS::Array{T}) where {T<:Real}

    # # Extract geometry from panels
    # NPAN = length(panels)
    # n_hats = zeros(NPAN,2)
    # t_hats = zeros(NPAN,2)
    # X = zeros(NPAN,2)
    # for n=1:NPAN
    #     n_hats[n,:] = panels[n].n_hat
    #     t_hats[n,:] = panels[n].t_hat
    #     X[n,1] = panels[n].R1[1]
    #     X[n,2] = panels[n].R1[2]
    # end

    X = [[point[1], point[2]] for point in boundary.bodyPTS]
    n_hats = [[n_hat[1], n_hat[2]] for n_hat in boundary.nHats]
    t_hats = [[t_hat[1], t_hat[2]] for t_hat in boundary.tHats]
    NPTS = boundary.NPTS_BODY

    # Guassian spreading and normalization constant
    sigma = 1.0

    # === Helper Function Definitions ===
    X_1i(i) = X[i,:]
    X_2i(i) = (i != NPTS) ? X[i+1,:] : X[1,:]
    Xj(j) = X[j,:]
    R_0i(i) = X_1i(i) - X_2i(i)
    R_1ij(i,j) = X_1i(i) - Xj(j)
    R_2ij(i,j) = X_2i(i) - Xj(j)
    r_0i(i) = norm(R_0i(i))
    r_1ij(i,j) = norm(R_1ij(i,j))
    r_2ij(i,j) = norm(R_2ij(i,j))
    function theta_ij(i,j)
        R1 = [R_1ij(i,j)[1],R_1ij(i,j)[2],0]
        R2 = [R_2ij(i,j)[1],R_2ij(i,j)[2],0]
        num = norm(cross(R1,R2))
        den = dot(R_1ij(i,j), R_2ij(i,j))
        return atan2(num,den)
    end

    function r_ij(i,j,thetaPrime)
        dotted = dot(R_1ij(i,j), t_hats[i,:])
        A = r_1ij(i,j)^3 * cos(thetaPrime)
        B = r_1ij(i,j) * cos(thetaPrime)*dotted^2
        C = sqrt(abs(r_1ij(i,j)^4 * dotted^2 * sin(thetaPrime)^2 - dotted^4 * sin(thetaPrime)^2))
        D = r_1ij(i,j)^2 * cos(thetaPrime)^2 - dotted^2
        return (A - B - C)/(D)
    end

    # === Coefficient Functions ===
    function Theta_ij(i,j)
        dottedT = dot(R_1ij(i,j), t_hats[i,:])
        if abs(dottedT) > 1e-14
            dottedN = dot(R_1ij(i,j), n_hats[i,:])
            b2 = dottedT^2 - r_1ij(i,j)^2
            aLim = r_1ij(i,j)
            bLim = sqrt(r_0i(i)^2 - 2*r_0i(i)*dottedT + r_1ij(i,j)^2)
            function f(r)
                f = (exp((-r^2)/(2*sigma^2))/r)*(1/sqrt(abs(b2 + r^2)))
                return f
            end
            return (dottedN/pi)*(quadgk(f,aLim,bLim)[1])
        else
            return 0
        end
    end

    function Lambda_ij(i,j)
        dotted = dot(R_1ij(i,j), t_hats[i,:])
        if abs(dotted) > 1e-14
            b2 = dotted^2 - r_1ij(i,j)^2
            aLim = r_1ij(i,j)
            bLim = sqrt(r_0i(i)^2 - 2*r_0i(i)*dotted + r_1ij(i,j)^2)
            g(r) = r*exp((-r^2)/(2*sigma^2))*(1/sqrt(b2 + r^2))
            return quadgk(g,aLim,bLim)[1]
        else
            return 0
        end
    end

    function phi_ij(i,j)
        R = norm(X_1i(j) - X_1i(i))
        return RBF_gauss(R)
        #return exp(-norm(X_1i(j) - X_1i(i))/(2*sigma^2))
    end

    # === Build coef matrix ===
    #matA = zeros(NPAN-1,NPAN-1)
    # for i=1:NPAN-1
    #     for j=1:NPAN-1
    #         if i == j
    #             matA[i,j] = phi_ij(i,j)
    #         else
    #             t = Theta_ij(i,j)
    #             l = Lambda_ij(i,j)
    #             matA[i,j] = phi_ij(i,j) - t + l
    #         end
    #     end
    # end

    #matA = [[(i == j) ? phi_ij(i,j) : phi_ij(i,j) - Theta_ij(i,j) + Lambda_ij(i,j) for j in 1:NPTS] for i in 1:NPTS]
    matA = zeros(NPTS, NPTS)
    for i=1:NPTS
        for j=1:NPTS
            if i == j
                matA[i,j] = phi_ij(i,j)
            else
                t = Theta_ij(i,j)
                l = Lambda_ij(i,j)
                matA[i,j] = phi_ij(i,j) - t + l
            end
        end
    end

    #return matA\RHS

    # === Solve Linear System ===
    alpha = matA\RHS
    return alpha



end

"""
`CalcVSCoefs(boundary, U_slip)`
"""
function CalcVSCoefs(boundary::Boundary; U_slip::Union{Nothing, Function}=nothing)

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
    rho_1(k) = CalcRho(boundary, X_k[k])                # ρ_1 of Eigendecomposition

    # === Coefficent functions ===
    function Theta_ki(k,i)
        res = 0
        for j in 1:NPTS
            if j != k
                res = res + (1/(2*pi)) * (dot(R_kj(k,j), n_hats[j]) * RBF_gauss(r_ji(j,i)) * del_S(j))/(r_kj(k,j)^2)
            end
        end
        return res
    end
    #Theta_ki(k,i) = sum([((k != j) ? (1/pi) * (dot(R_kj(k,j), n_hats[j]) * RBF_gauss(r_ji(j,i)) * del_S(j))/(r_kj(k,j)^2) : 0) for j in 1:NPTS] )
    Lambda_ki(k,i) = sum( [RBF_gauss(r_ji(j,i)) * (rho_1(k))/(L) for j in 1:NPTS] )
    phi_ki(k,i) = RBF_gauss(r_ki(k,i))

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
                X_eval::Array{T}) where {T<:Real}

    # === Unpack geometry ===
    XP = [[point[1], point[2], 0.0] for point in boundary.bodyPTS]
    NPTS = boundary.NPTS_BODY

    # === Summation over all points for gamma ===
    gamma = sum( [alpha[i] * RBF_gauss(norm(X_eval - XP[i])) for i in 1:NPTS ] )
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





    # # === Integrand function ===
    # function F(theta)
    #     rP_hat = [cos(theta), sin(theta), 0]
    #     X_P = radius .* rP_hat
    #     gamma = CalcVS(boundary, alpha, X_P) .* [0, 0, 1]
    #     num = cross(X_eval - X_P, gamma) * radius
    #     den = 4 * pi * norm(X_eval - X_P)^3
    #     return num ./ den
    # end
    #
    # Fx(theta) = F(theta)[1]
    # Fy(theta) = F(theta)[2]
    # Fz(theta) = F(theta)[3]
    #
    # # === Evaluate integral ===
    # rtol = 1e-3
    # maxevals = 1e3
    #
    # return [quadgk(Fx, 0, 2*pi, rtol=rtol, maxevals=maxevals)[1],
    #         quadgk(Fy, 0, 2*pi, rtol=rtol, maxevals=maxevals)[1],
    #         quadgk(Fz, 0, 2*pi, rtol=rtol, maxevals=maxevals)[1]]
end

"""
`CalcVortexSheetVelocity(panels,alpha,X,sigma=0.2,a=1.0)`

Given a geometry of panels, the strength of the vortex particle on the panels,
and an X location, calculates the velocity induced by the vortex sheet at that point.

# ARGUMENTS
* `panels::Array{Panel2D}` : Array of Panel2D objects representing a body
* `alpha::Array{Float64`   : Radial basis function coefficients previously solved for
* `X::Array{Float64}`      : Point to compute velocity
* `sigma::Float64=0.2`     : Guassian width
* `a::Float64=1.0`         : Normalization constant for Guassian

# OUTPUTS
* `u_vor::Float64`         : Calculated velocity from vortex sheet
"""
function CalcVortexSheetVelocity(panels::Array{SciTools.LineSegment},
                                 alpha::Array{Float64},
                                 X::Array{Float64};
                                 smoothRadius=1e-6,
                                 sigma::Float64=0.2,
                                 a::Float64=1.0)

    # Extract geometry from panels
    NPAN = length(panels)
    XP = zeros(NPAN,3)
    for n=1:NPAN
        XP[n,1] = panels[n].R1[1]
        XP[n,2] = panels[n].R1[2]
    end

    # Geometry function definitions
    R0(X1,X2) = X1 - X2
    R1(X,X1) = X - X1
    R2(X,X2) = X - X2
    r0(X1,X2) = norm(R0(X1,X2))
    r1(X,X1) = norm(R1(X,X1))
    r2(X,X2) = norm(R1(X,X2))

    # Velocity function
    function U(X,X1,X2,Gamma)
        A = Gamma/(4*pi)
        B = cross(R1(X,X1),R2(X,X2))/(norm(cross(R1(X,X1),R2(X,X2)))^2)
        C = (dot(R0(X1,X2),R1(X,X1)))/(r1(X,X1))
        D = (dot(R0(X1,X2),R2(X,X2)))/(r2(X,X2))
        u = A * B * (C - D)
        return u
    end

    U_VS = zeros(3,1)
    for i=1:NPAN-1
        gamma = CalcVortexSheet(panels,alpha,X,sigma=sigma,a=a)
        toAdd = U(X,XP[i,:],XP[i+1,:],gamma*panels[i].L)
        if isnan.(toAdd[1]) == false
            U_VS = U_VS + toAdd
        end

    end

    return U_VS
end
