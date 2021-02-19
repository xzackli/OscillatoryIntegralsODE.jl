using OscillatoryIntegralsODE
using SpecialFunctions
using BenchmarkTools
using PyPlot
using OrdinaryDiffEq

nu = 100.
r = 100.
f(x) = exp(-x^2/16)

clf()
xs = 0:0.001:5.0
plot(xs, [f(x) * besselj(nu, r * x) for x in xs], "-", label=raw"$e^{-x^2/16} J_{100}(100x)$")
legend()
gcf()
# savefig("docs/src/assets/bessel.svg")

##
bi = BesselJIntegral{Float64}(nu, r)
result = levintegrate(bi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)

ref = 0.006311630277650583
print("relative error = ", (result - ref) / ref, "\n")

##
@btime levintegrate($bi, $f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
