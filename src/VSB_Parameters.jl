"""

"""
mutable struct Parameters

    boundary::Boundary
    U_field::Union{Nothing, Function}
    alphas::Array{T} where {T<:Real}
    etas::Array{T} where {T<:Real}
    verbose::Bool

    function Parameters(boundary, U_field; verbose = false)

        # === Calculate etas ===
        if verbose
            println("   CALCULATING ETAS IN PARAMETERS...")
        end
        etas = CalcRhoCoefs(boundary)

        # === Calculate alphas ===
        if verbose
            println("   CALCULATING ALPHAS IN PARAMETERS...")
        end
        alphas = CalcVSCoefs(boundary, etas, U_slip = U_field)

        new(boundary, U_field, alphas, etas)

    end

end
