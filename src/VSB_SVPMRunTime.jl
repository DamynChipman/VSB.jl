"""
`SVPMRunTime(pfield, t, dt)`

[description]

Psuedocode:
SVPMRunTime
    Add particles on body
    Calc: U_slip = U(pfield) + U_inf
    Calc VS Strength
    Diffuse vorticity
end

# ARGUMENTS
* `pfield::ParticleField`   : SimpleVPM ParticleField
* `t::Real`                 : Current time
* `dt::Real`                : Time step
* `body::Array{Array{T,1}}` : Matrix of points representing physical boundary

# OUTPUTS
* `breakflag::Bool`       : Flag for stopping SVPM simulation
"""
function SVPMRunTime(pfield::SVPM.ParticleField,
                     t::Real,
                     dt::Real,
                     body::Array{Array{T,1},1}) where{T<:Real}

    # === Add Particles to field around Boundary ===



    # === Calculate U_slip ===



    # === Calculate Vortex Sheet strength ===



    # === Diffuse vorticity ===



    # === Returns ===

    
    return true
end
