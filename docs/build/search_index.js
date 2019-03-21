var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "VSB.jl Documentation",
    "title": "VSB.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "#VSB.jl-Documentation-1",
    "page": "VSB.jl Documentation",
    "title": "VSB.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "#Vortex-Sheet-Boundary-Method-for-Fluid-Systems-1",
    "page": "VSB.jl Documentation",
    "title": "Vortex Sheet Boundary Method for Fluid Systems",
    "category": "section",
    "text": "Author: Damyn ChipmanA mesh-less and Lagrangian approach to address the presence of boundaries for application with the Vortex Particle Method."
},

{
    "location": "#Install-1",
    "page": "VSB.jl Documentation",
    "title": "Install",
    "category": "section",
    "text": "pkg> add https://github.com/camperD/VSB.jl"
},

{
    "location": "#Orientation-1",
    "page": "VSB.jl Documentation",
    "title": "Orientation",
    "category": "section",
    "text": "Theory for the method, including a detailed mathematical derivation and numerical discretization, can be found under TheoryFor usage examples, see [TODO: Finish and add tests. Link to here]."
},

{
    "location": "guide/#",
    "page": "Usage Guide",
    "title": "Usage Guide",
    "category": "page",
    "text": ""
},

{
    "location": "guide/#Usage-Guide-1",
    "page": "Usage Guide",
    "title": "Usage Guide",
    "category": "section",
    "text": "This is a guide going through some of the test cases for Vortex Sheet Boundary.using VSB         # Import the package\n\nvar = 3"
},

{
    "location": "docstrings/#",
    "page": "Code Documentation",
    "title": "Code Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "docstrings/#VSB.Boundary",
    "page": "Code Documentation",
    "title": "VSB.Boundary",
    "category": "type",
    "text": "Boundary(body_pts,t_hats,n_hats)\n\nStruct containing boundary information. The user must supply the location of the body points (body_pts), and the tangent and normal vectors of the body at each body point.\n\nARGUMENTS\n\nbody_pts   : Array of the form [[x1,y1,z1], [x2,y2,z2], ...] containing the                locations of the boundary nodes. Assumes a CCW orientation.\nt_hats     : Array of the same form containing the tangent unit vectors at                each ponit of the boundary. Same size as body_pts.\nn_hats     : Array of the same form containing the normal unit vectors at                each ponit of the boundary. Same size as body_pts.\n\nPROPERTIES\n\nbodyPTS    : Array of the form [[x1,y1,z1], [x2,y2,z2], ...] containing the                locations of the boundary nodes.\ntHats      : Array of the same form containing the tangent unit vectors at                each ponit of the boundary. Same size as body_pts.\nnHats      : Array of the same form containing the normal unit vectors at                each ponit of the boundary. Same size as body_pts.\nNPTS_BODY  : Number of points in the body.\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.Parameters",
    "page": "Code Documentation",
    "title": "VSB.Parameters",
    "category": "type",
    "text": "Parameters(boundary, U_field)\n\nStruct that represents problem parameters that are determined by the geometry. A different Parameters object should be created (or a single one updated) for every geometry (boundary) or whenever the velocity field (U_field) changes.\n\nCalls CalcRhoCoefs and CalcVSCoefs to update the parameters.\n\nARGUMENTS\n\nboundary    : Boundary object representing geometry of boundary.\nU_field     : Function of the form U_field(X). Gives the velocity at every                 point in the domain for point X.\nsigma       : Gaussian width. DEFAULT = 0.2.\nverbose     : Boolean for printing status of calculation.\n\nPROPERTIES\n\nboundary    : Boundary object representing geometry of boundary.\nU_field     : Function of the form U_field(X). Gives the velocity at every                 point in the domain for point X.\nalphas      : Array containing RBF coefficients for vortex sheet strength.\netas        : Array containing RBF coefficients for geometric eigen decomposition                 function.\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.CalcRhoCoefs",
    "page": "Code Documentation",
    "title": "VSB.CalcRhoCoefs",
    "category": "function",
    "text": "`CalcRhoCoefs(self::Boundary)`\n\nTakes a Boundary object and calculates the RBF coefficients for the eigenfunction of the kernel decomposition.\n\nARGUMENTS\n\nself::Boundary   : Boundary object\nsigma = 0.2      : Gaussian width\n\nRETURNS\n\netas::Array{Float64} : Array containing the RBF coefficient values for rho\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.CalcRho",
    "page": "Code Documentation",
    "title": "VSB.CalcRho",
    "category": "function",
    "text": "`CalcRho(self::Boundary, etas, X_eval)`\n\nEvaluates the RBF approxiamation function: rho(X) = sum(etai * phi(X - X)).\n\nARGUMENTS\n\nself::Boundary          : Boundary object\netas::Array{Float64}    : Array containing RBF coefficients\nX_eval::Array{Float64}  : Point to evaluate function at\nsigma = 0.2             : Gaussian width\n\nRETURNS\n\nrho::Float64            : Evaluated function value\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.Circle",
    "page": "Code Documentation",
    "title": "VSB.Circle",
    "category": "function",
    "text": "`Circle(R, O, N)`\n\nCreates a Boundary object in the shape of a circle of radius R at origin O.\n\nARGUMENTS\n\nR::Float64        : Radius of cirlce\nO::Array{Float64} : Origin of cirlce of form [X, Y, Z]\nN::Int64          : Number of points\n\nRETURNS\n\ncircle::Boundary  : Boundary object of cirlce\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.NACA4",
    "page": "Code Documentation",
    "title": "VSB.NACA4",
    "category": "function",
    "text": "`NACA4(numb, N, c)`\n\nNACA Four-Digit Series Airfoil. Returns two lists of ordered pairs representing a NACA Four-Digit Airfoil.\n\nARGUMENTS\n\nnumb::String     : Four digit series. Symmetric airfoil given by \"00xx\"\nN::Int64         : Number of points\nc::Float64=1.0   : Length of airfoil. Defaults to 1.0 for appropiate scaling\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.RBF_gauss",
    "page": "Code Documentation",
    "title": "VSB.RBF_gauss",
    "category": "function",
    "text": "`RBF_gauss(R;A=1.0,sigma=0.2,deriv=0)`\n\nKernel function for Gaussian RBF Interpolation. Returns the magnitude of evaluated function.\n\nThe optional argument deriv is for the order of the derivative desired, meaning 0 for base, 1 for first derivative (gradient) and 2 for second derivative (Laplacian).\n\nARGUMENTS\n\nR::Float64 : Magnitude of radius vector -> R = |Xj - Xi|\nA=1.0      : Normalization value\nsigma=0.2  : Gaussian spreading\nderiv=0    : Order of desired derivative\n\nOUTPUTS\n\nres        : Calculated magnitude\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.CalcVSCoefs",
    "page": "Code Documentation",
    "title": "VSB.CalcVSCoefs",
    "category": "function",
    "text": "`CalcVSCoefs(boundary, etas; U_slip, sigma)`\n\nCalculates the vortex sheet RBF coefficients.\n\nARGUMENTS\n\nboundary::Boundary                       : Boundary object\netas::Array{T}                           : RBF coefficients for rho\nU_slip::Union{Nothing, Function}=nothing : Slip velocity field function: U_slip(X)\nsigma = 0.2                              : Gaussian width\n\nwhere {T<:Real}\n\nRETURNS\n\nalphas::Array{Float64}                   : Array of RBF coefficients\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.CalcVS",
    "page": "Code Documentation",
    "title": "VSB.CalcVS",
    "category": "function",
    "text": "`CalcVS(boundary, alpha, X_eval)`\n\nEvaluates the RBF interpolation function: gamma(X) = sum(alphai * phi(X - Xi)).\n\nARGUMENTS\n\nboundary::Boundary : Boundary object\nalpha::Array{T}    : Array containing RBF coefficients\nX_eval::Array{T}   : Point to evaluate function at\nsigma = 0.2        : Gaussian width\n\nwhere {T<:Real}\n\nRETURNS\n\ngamma::Float64     : Evaluated function value\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.CalcVSDiffusionCoefs",
    "page": "Code Documentation",
    "title": "VSB.CalcVSDiffusionCoefs",
    "category": "function",
    "text": "`CalcVSDiffusionCoefs(boundary, alphas, dt, nu)`\n\nCalculates the RBF coefficients for the diffusion step.\n\nARGUMENTS\n\nboundary::Boundary     : Boundary object\netas::Array{T}         : RBF coefficients for rho\ndt::Real               : Time step\nnu::Real               : Kinematic viscosity\nsigma = 0.2            : Gaussian width\n\nwhere {T<:Real}\n\nRETURNS\n\nbetas::Array{Float64}  : Array of RBF coefficients\n\n\n\n\n\n"
},

{
    "location": "docstrings/#VSB.CalcDiffusion",
    "page": "Code Documentation",
    "title": "VSB.CalcDiffusion",
    "category": "function",
    "text": "`CalcDiffusion(boundary, beta, X_eval)`\n\nEvaluates the RBF interpolation function: omega(X) = sum(betai * phi(X - Xi)).\n\nARGUMENTS\n\nboundary::Boundary : Boundary object\nbeta::Array{T}     : Array containing RBF coefficients\nX_eval::Array{T}   : Point to evaluate function at\nsigma = 0.2        : Gaussian width\n\nwhere {T<:Real}\n\nRETURNS\n\nomega::Float64     : Evaluated function value\n\n\n\n\n\n"
},

{
    "location": "docstrings/#Code-Documentation-1",
    "page": "Code Documentation",
    "title": "Code Documentation",
    "category": "section",
    "text": "Below are all the available functions and structs with their associated docstring detailing their usage.Boundary\nParameters\nCalcRhoCoefs\nCalcRho\nCircle\nNACA4\nRBF_gauss\nCalcVSCoefs\nCalcVS\nCalcVSDiffusionCoefs\nCalcDiffusion"
},

]}
