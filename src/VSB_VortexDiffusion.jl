"""
`VortexDiffusion()`

[description]

# ARGUMENTS
* () : []

# OUTPUTS
* () : []

"""
function VortexDiffusion(pfield::SVPM.ParticleField,
                         boundary::Boundary,
                         X_eval::Array{T},
                         dt::Real,
                         RHS::Array{T};
                         theta_CN::Real=0.5) where {T<:Real}


    # Constants
    N_BODY = boundary.NPTS_BODY
    N_FIELD = length(pfield.particles)
    CONST1 = pfield.nu * dt * theta_CN
    X = [[point[1], point[2]] for point in boundary.bodyPTS]

    # Helper functions
    rHat(X1,X2) = (X1 - X2)/(norm(X1 - X2))
    r(X1,X2) = norm(X1 - X2)

    # Build matrix
    A = zeros(N_FIELD,N_FIELD)

    # NxM Matrix (upper half)
    for i = 1:N_BODY
        for j = 1:N_FIELD
            Xi = pfield.particles[i].X
            Xj = pfield.particles[j].X
            R = r(Xi,Xj)
            RHAT = rHat(Xi,Xj)
            A[i,j] = dot(RBF_gauss(R,deriv=1).*RHAT, boundary.nHats[i])
            println("A[",i,", ",j,"] = ",A[i,j])
        end
    end

    # (M-N+1)xM matrix (lower half)
    for i = (N_BODY+1):N_FIELD
        for j = 1:N_FIELD
            Xi = pfield.particles[i].X
            Xj = pfield.particles[j].X
            R = r(Xi,Xj)
            A[i,j] = RBF_gauss(R,deriv=0) - CONST1 * RBF_gauss(R,deriv=2)
            println("A[",i,", ",j,"] = ",A[i,j])
        end
    end

    # Solve system for coefs
    beta = A\RHS
    println("BETA INSIDE DIFFUSION:")
    println(beta)

    omega = 0
    for i=1:N_BODY
        omega = omega + beta[i] * RBF_gauss(norm(X_eval[1:2] - X[i]))
    end
    return omega

end
