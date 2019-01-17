"""
======================== VORTEX SHEET BOUNDARY =================================
HEADER FILE

@author: Damyn Chipman


"""
module VSB

# ===== Imports =====
using SciTools
using QuadGK

# ===== Files =====
file_names = ["CalcVortexSheet", "Verification1"]
for header_name in file_names
    include("VSB_"*header_name*".jl")
end

# ========================= END OF VSB.jl MODULE ===============================
end
