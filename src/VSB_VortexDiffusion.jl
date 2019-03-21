"""
`VortexDiffusion()`

[description]

# ARGUMENTS
* () : []

# OUTPUTS
* () : []

"""
function CalcDiffusionCoefs(pfield::SVPM.ParticleField,
                            boundary::Boundary,
                            alphas::Array{T},
                            dt::Real;
                            theta_CN::Real=0.5) where {T<:Real}


    # === Extract geometry from boundary ===
    X = [[point[1], point[2], 0.0] for point in boundary.bodyPTS]
    n_hats = [[n_hat[1], n_hat[2]] for n_hat in boundary.nHats]
    t_hats = [[t_hat[1], t_hat[2]] for t_hat in boundary.tHats]
    NBODY = boundary.NPTS_BODY

    # === Constants ===
    CONST1 = pfield.nu * dt * theta_CN
    CONST2 = pfield.nu * dt

    # === Set of points: Union of points from surface and particle field ===
    XPTS = []
    for point in X
        push!(XPTS, point)
    end
    if length(pfield.particles) != 0
        for particle in pfield.particles
            push!(XPTS, particle.X)
        end
    end
    NFIELD = length(XPTS)

    # === Build matrix ===
    matA = zeros(NFIELD, NFIELD)

    # NxM Matrix (upper half)
    for i = 1:NBODY
        for j = 1:NFIELD
            if i != j
                R = XPTS[j] - XPTS[i]
                r = norm(R)
                r_hat = R ./ r
                matA[i,j] = dot(RBF_gauss(r, deriv=1) .* r_hat, boundary.nHats[i])
                #println(" Xi = ",Xi," ... Xj = ",Xj," ... R = ",R," ... A[",i,", ",j,"] = ",A[i,j])
            else
                R = [0.0, 0.0, 0.0]
                r = 0.0
                r_hat = [0.0, 0.0, 0.0]
                matA[i,j] = dot(RBF_gauss(r, deriv=1) .* r_hat, boundary.nHats[i])
                #println("I == J ... A[",i,", ",j,"] = ",A[i,j])
            end

        end
    end

    # (M-N+1)xM matrix (lower half)
    for i = (NBODY+1):NFIELD
        for j = 1:NFIELD
            if i != j
                R = XPTS[j] - XPTS[i]
                r = norm(R)
                matA[i,j] = RBF_gauss(r, deriv=0) - CONST1 * RBF_gauss(r, deriv=2)
                #println(" Xi = ",Xi," ... Xj = ",Xj," ... R = ",R," ... A[",i,", ",j,"] = ",A[i,j])
            else
                R = [0.0, 0.0, 0.0]
                r = 0.0
                matA[i,j] = RBF_gauss(r, deriv=0) - CONST1 * RBF_gauss(r, deriv=2)
                #println("I == J ... A[",i,", ",j,"] = ",A[i,j])
            end
        end
    end

    # === Build RHS ===
    RHS = zeros(NFIELD)
    for i = 1:NBODY
        RHS[i] = CalcVS(boundary, alphas, X[i]) / CONST2
    end

    # === Solve system for coefs ===
    beta = matA\RHS
    return beta

end

"""
    `CalcVSDiffusionCoefs(boundary, alphas, dt, nu)`

Calculates the RBF coefficients for the diffusion step.

# ARGUMENTS
* `boundary::Boundary`     : Boundary object
* `etas::Array{T}`         : RBF coefficients for rho
* `dt::Real`               : Time step
* `nu::Real`               : Kinematic viscosity
* `sigma = 0.2`            : Gaussian width
where {T<:Real}

# RETURNS
* `betas::Array{Float64}`  : Array of RBF coefficients
"""
function CalcVSDiffusionCoefs(boundary::Boundary,
                              alphas::Array{T},
                              dt::Real,
                              nu::Real;
                              sigma=0.2) where {T<:Real}

    # === Extract geometry from boundary ===
    X = [[point[1], point[2], 0.0] for point in boundary.bodyPTS]
    n_hats = [[n_hat[1], n_hat[2], 0.0] for n_hat in boundary.nHats]
    t_hats = [[t_hat[1], t_hat[2], 0.0] for t_hat in boundary.tHats]
    NPTS = boundary.NPTS_BODY

    # === Constants ===
    CONST1 = -1 / (nu * dt)

    # === Build matrix ===
    A = zeros(NPTS, NPTS)
    for k in 1:NPTS
        for i in 1:NPTS
            if k != i
                R = X[k] - X[i]
                r = norm(R)
                r_hat = R ./ r
                n_hat = n_hats[k]
                A[k,i] = dot(RBF_gauss(r, deriv = 1, sigma=sigma) .* r_hat, n_hat)
            else
                R = [0, 0, 0]
                r = norm(R)
                r_hat = [0, 0, 0]
                n_hat = n_hats[k]
                A[k,i] = dot(RBF_gauss(r, deriv = 1, sigma=sigma) .* r_hat, n_hat)
            end
        end
    end

    # === Build RHS ===
    b = zeros(NPTS)
    for k in 1:NPTS
        b[k] = CalcVS(boundary, alphas, X[k]) * CONST1
    end

    # === Solve linear system ===
    betas = A\b
    return betas

end


"""
    `CalcDiffusion(boundary, beta, X_eval)`

Evaluates the RBF interpolation function: omega(X) = sum(beta_i * phi(X - X_i)).

# ARGUMENTS
* `boundary::Boundary` : Boundary object
* `beta::Array{T}`     : Array containing RBF coefficients
* `X_eval::Array{T}`   : Point to evaluate function at
* `sigma = 0.2`        : Gaussian width
where {T<:Real}

# RETURNS
* `omega::Float64`     : Evaluated function value
"""
function CalcDiffusion(boundary::Boundary,
                       beta::Array{T},
                       X_eval::Array{T};
                       sigma=0.2) where {T<:Real}

    # === Unpack geometry ===
    XP = [[point[1], point[2], 0.0] for point in boundary.bodyPTS]
    NPTS = boundary.NPTS_BODY

    # === Summation over all points for gamma ===
    omega = sum( [beta[i] * RBF_gauss(norm(X_eval - XP[i]), sigma=sigma) for i in 1:NPTS ] )
    return omega
end
