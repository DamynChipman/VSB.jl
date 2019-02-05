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
                         dt::Real,
                         RHS::Array{T,1};
                         theta_CN::Real=0.5) where {T<:Real}


    # Constants
    N_BODY = boundary.NPTS_BODY
    N_FIELD = length(pfield.particles)
    CONST1 = pfield.nu * dt * theta_CN

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
        end
    end

    # (M-N+1)xM matrix (lower half)
    for i = (N_BODY+1):N_FIELD
        for j = 1:N_FIELD
            Xi = pfield.particles[i].X
            Xj = pfield.particles[j].X
            R = r(Xi,Xj)
            A[i,j] = RBF_guass(R,deriv=0) - CONST1 * RBF_guass(R,deriv=2)
        end
    end

    # Solve system for coefs
    coefs = A\RHS
    return coefs
end
