"""
    Parameters(boundary, U_field)

Struct that represents problem parameters that are determined by the geometry. A
different Parameters object should be created (or a single one updated) for every
geometry (boundary) or whenever the velocity field (U_field) changes.

Calls CalcRhoCoefs and CalcVSCoefs to update the parameters.

# ARGUMENTS
* `boundary`    : Boundary object representing geometry of boundary.
* `U_field`     : Function of the form U_field(X). Gives the velocity at every
                  point in the domain for point X.
* `sigma`       : Gaussian width. DEFAULT = 0.2.
* `verbose`     : Boolean for printing status of calculation.

# PROPERTIES
* `boundary`    : Boundary object representing geometry of boundary.
* `U_field`     : Function of the form U_field(X). Gives the velocity at every
                  point in the domain for point X.
* `alphas`      : Array containing RBF coefficients for vortex sheet strength.
* `etas`        : Array containing RBF coefficients for geometric eigen decomposition
                  function.
"""
mutable struct Parameters

    boundary::Boundary
    U_field::Union{Nothing, Function}
    alphas::Array{T} where {T<:Real}
    etas::Array{T} where {T<:Real}
    verbose::Bool

    function Parameters(boundary, U_field; sigma = 0.2, verbose = false)

        # === Calculate etas ===
        if verbose
            println("   CALCULATING ETAS IN PARAMETERS...")
        end
        etas = CalcRhoCoefs(boundary, sigma=sigma)

        # === Calculate alphas ===
        if verbose
            println("   CALCULATING ALPHAS IN PARAMETERS...")
        end
        alphas = CalcVSCoefs(boundary, etas, U_slip = U_field, sigma=sigma)

        new(boundary, U_field, alphas, etas)

    end

end
