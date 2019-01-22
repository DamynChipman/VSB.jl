"""
======================== VORTEX SHEET BOUNDARY =================================
HEADER FILE

Vortex Sheet Boundary method for coupling with Vortex Particle Method. Written
Julia 1.0

@author: Damyn Chipman

"""
module VSB

# ===== Imports =====
using SciTools
using QuadGK
using LinearAlgebra
using SimpleVPM
SVPM = SimpleVPM

# ===== GLOBAL VARIABLES =====
const NPTS_BODY                             # Number of poitns on vortex body
const module_path = splitdir(@__file__)[1]  # File path to module VSB.jl

# ===== Files =====
file_names = ["CalcVortexSheet","RBF", "VortexDiffusion", "SVPMRunTime"]
for header_name in file_names
    include("VSB_"*header_name*".jl")
end

# ========================= END OF VSB.jl MODULE ===============================
end
