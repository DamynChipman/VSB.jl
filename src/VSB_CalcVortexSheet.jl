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
function CalcVortexSheetCoef(panels::Array{LineSegment},RHS::Array{Float64})

    # Extract geometry from panels
    NPAN = length(panels)
    n_hats = zeros(NPAN,2)
    t_hats = zeros(NPAN,2)
    X = zeros(NPAN,2)
    for n=1:NPAN
        n_hats[n,:] = panels[n].n_hat
        t_hats[n,:] = panels[n].t_hat
        X[n,1] = panels[n].r1[1]
        X[n,2] = panels[n].r1[2]
    end

    # Guassian spreading and normalization constant
    a = 1
    sigma = 0.2

    # Function definitions
    X_1i(i) = X[i,:]
    X_2i(i) = X[i+1,:]
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
        den = dot(2,R_1ij(i,j),1,R_2ij(i,j),1)
        return atan2(num,den)
    end

    function r_ij(i,j,thetaPrime)
        dotted = dot(2,R_1ij(i,j),1,t_hats[i,:],1)
        A = r_1ij(i,j)^3 * cos(thetaPrime)
        B = r_1ij(i,j) * cos(thetaPrime)*dotted^2
        C = sqrt(abs(r_1ij(i,j)^4 * dotted^2 * sin(thetaPrime)^2 - dotted^4 * sin(thetaPrime)^2))
        D = r_1ij(i,j)^2 * cos(thetaPrime)^2 - dotted^2
        return (A - B - C)/(D)
    end

    # --- Coefficient Functions ---
    function Theta_ij(i,j)
        dottedT = dot(2,R_1ij(i,j),1,t_hats[i,:],1)
        if abs(dottedT) > 1e-14
            dottedN = dot(2,R_1ij(i,j),1,n_hats[i,:],1)
            b2 = dottedT^2 - r_1ij(i,j)^2
            aLim = r_1ij(i,j)
            bLim = sqrt(r_0i(i)^2 - 2*r_0i(i)*dottedT + r_1ij(i,j)^2)
            function f(r)
                f = (exp((-r^2)/(2*sigma^2))/r)*(1/sqrt(abs(b2 + r^2)))
                return f
            end
            return (a*dottedN/pi)*(quadgk(f,aLim,bLim)[1])
        else
            return 0
        end
    end

    function Lambda_ij(i,j)
        dotted = dot(2,R_1ij(i,j),1,t_hats[i,:],1)
        if abs(dotted) > 1e-14
            b2 = dotted^2 - r_1ij(i,j)^2
            aLim = r_1ij(i,j)
            bLim = sqrt(r_0i(i)^2 - 2*r_0i(i)*dotted + r_1ij(i,j)^2)
            g(r) = a*r*exp((-r^2)/(2*sigma^2))*(1/sqrt(b2 + r^2))
            return quadgk(g,aLim,bLim)[1]
        else
            return 0
        end
    end

    function phi_ij(i,j)
        return a*exp(-norm(X_1i(j) - X_1i(i))/(2*sigma^2))
    end

    # Build coef matrix
    matA = zeros(NPAN-1,NPAN-1)
    for i=1:NPAN-1
        for j=1:NPAN-1
            if i == j
                matA[i,j] = phi_ij(i,j)
            else
                t = Theta_ij(i,j)
                l = Lambda_ij(i,j)
                matA[i,j] = phi_ij(i,j) - t + l
            end
        end
    end

    return matA\RHS
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
function CalcVortexSheet(panels::Array{LineSegment},
                         alpha::Array{Float64},
                         X::Array{Float64};
                         sigma::Float64=0.2,
                         a::Float64=1.0)
    gamma = 0
    N = length(alpha)
    for i=1:N
        gamma = gamma + a*alpha[i]*exp((-norm(X[1:2] - panels[i].r1)^2)/(2*sigma^2))
    end
    return gamma
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
function CalcVortexSheetVelocity(panels::Array{LineSegment},
                                 alpha::Array{Float64},
                                 X::Array{Float64};
                                 smoothRadius=1e-6,
                                 sigma::Float64=0.2,
                                 a::Float64=1.0)

    # Extract geometry from panels
    NPAN = length(panels)
    XP = zeros(NPAN,3)
    for n=1:NPAN
        XP[n,1] = panels[n].r1[1]
        XP[n,2] = panels[n].r1[2]
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
        C = (dot(3,R0(X1,X2),1,R1(X,X1),1))/(r1(X,X1))
        D = (dot(3,R0(X1,X2),1,R2(X,X2),1))/(r2(X,X2))
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
