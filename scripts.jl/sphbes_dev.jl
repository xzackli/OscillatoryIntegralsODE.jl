using OscillatoryIntegralsODE
using SpecialFunctions
using StaticArrays
using BenchmarkTools
using PyPlot
using OrdinaryDiffEq

nu = 100
r = 100
f(x) = exp(-x^2/16)

clf()
xs = 0:0.001:5.0
plot(xs, [f(x) * sphericalbesselj(nu, r * x) for x in xs], "-", label=raw"$e^{-x^2/16} j_{100}(100x)$")
legend()
gcf()

##
si = SphericalBesselJIntegral{Float64}(nu, r)
result = levintegrate(si, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)

ref = 0.0008322179291456167
print("relative error = ", (result - ref) / ref, "\n")

## timing
@btime levintegrate($si, $f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
