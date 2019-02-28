# ===== VSB Verifications: RBF Gaussian Plots =====

println("=== BEGINNING GaussPlotting.jl ===")

# == Imports ==
using VSB
using Plots
pyplot()
using LaTeXStrings

cd("/Users/Damyn/Documents/BYU/FLOW Lab/VSB/verifications")

println("   PLOTTING...")
Xg = -1:0.01:1
Yg1 = [VSB.RBF_gauss(X,A=1.0,deriv = 0,sigma=0.2) for X in Xg]
Yg2 = [VSB.RBF_gauss(X,A=0.33 ,deriv = 1,sigma=0.2) for X in Xg]
Yg3 = [VSB.RBF_gauss(X,A=0.1153,deriv = 2,sigma=0.2) for X in Xg]
plt = plot(Xg, Yg1, line = (:blue), label = L"\phi(x)")
plot!(Xg,Yg2, line = (:red), label = L"\nabla \phi(x)")
plot!(Xg,Yg3, line = (:green), label = L"∇^2 ϕ(x)")
savefig(plt,"VSB_Verification_RBFGauss.pdf")

println("=== END OF GaussPlotting.jl ===")
# L"\nabla^2 \phi(|\textbf{x} - \textbf{x}_i|), A = 0.5"
