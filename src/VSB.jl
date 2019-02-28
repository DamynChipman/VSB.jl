"""
======================== VORTEX SHEET BOUNDARY =================================
HEADER FILE

Vortex Sheet Boundary method for coupling with Vortex Particle Method. Written
Julia 1.0

See README.md for details on how to install and docs/Theory.pdf for derivation.

@author: Damyn Chipman - 2019
================================================================================
"""
module VSB

# ===== Imports =====
using SciTools
using QuadGK
using LinearAlgebra
using SimpleVPM
SVPM = SimpleVPM

# ===== GLOBAL VARIABLES =====
const NUMB_MATRIX = Array{Array{T,1},1} where {T<:Real} # Numerical matrix type
#const module_path = splitdir(@__file__)[1]             # File path to module VSB.jl

# ===== Files =====
file_names = ["RBF",
              "SVPMRunTime",
              "Boundary",
              "CalcVortexSheet",
              "VortexDiffusion",
              "Parameters"]
for header_name in file_names
    include("VSB_"*header_name*".jl")
end

# ========================= END OF VSB.jl MODULE ===============================
end
