"""

"""
mutable struct Parameters

    boundary::Boundary
    U_field::Union{Nothing, Function}
    alphas::Array{T} where {T<:Real}
    etas::Array{T} where {T<:Real}

    function Parameters(boundary, U_field)

        etas = CalcRhoCoefs(boundary)
        alphas = CalcVSCoefs(boundary, etas, U_slip = U_field)

        new(boundary, U_field, alphas, etas)

    end

end
