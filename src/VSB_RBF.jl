"""
`RBF_gauss(R;A=1.0,sigma=1.0,deriv=0)`

Kernel function for Gaussian RBF Interpolation. Returns the magnitude of evaluated
function.

The optional argument `deriv` is for the order of the derivative desired, meaning
0 for base, 1 for first derivative (gradient) and 2 for second derivative (Laplacian)

# ARGUMENTS
* `R::Float64` : Magnitude of radius vector -> R = |X_j - X_i|
* `A=1.0`      : Normalization value
* `sigma=1.0`  : Gaussian spreading
* `deriv=0`    : Order of desired derivative

# OUTPUTS
* `res`        : Calculated magnitude
"""
function RBF_gauss(R::Float64;
                   A=1.0,sigma=0.2,deriv=0)

    if deriv == 0
        const const1 = 1/(2*pi)^(3/2)
        res = const1*exp(-R^2/2)

        #res = SVPM.zeta_gauserf(R/sigma)
        #res = A*exp(-(R^2)/(2*sigma^2))
    elseif deriv == 1
        res = -((A*R)/(sigma^2))*exp(-(R^2)/(2*sigma^2))
    elseif deriv == 2
        res = ((A/sigma)^2)*((R/sigma)^2 - 3)*exp(-(R^2)/(2*sigma^2))
    end

    return res
end
