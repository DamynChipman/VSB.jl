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
* `pfield::ParticleField` : SimpleVPM ParticleField
* `t::Real`               : Current time
* `dt::Real`              : Time step

# OUTPUTS
* `breakflag::Bool`       : Flag for stopping SVPM simulation
"""
function SVPMRunTime(pfield::SVPM.ParticleField, t::Real, dt::Real)

end
