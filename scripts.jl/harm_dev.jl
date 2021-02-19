using OscillatoryIntegralsODE
using SpecialFunctions
using StaticArrays
using BenchmarkTools
using PyPlot
using OrdinaryDiffEq

ω = 100.
f(x) = exp(-x^2/16)

clf()
xs = 0:0.001:5.0
plot(xs, real.([f(x) * exp(1im * ω * x) for x in xs]), "-", label=raw"$e^{-x^2/16} e^{-i 100 x}$")
legend()
gcf()

##
hi = HarmonicIntegral{Complex{Float64}}(ω)
result = levintegrate(hi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)

ref = 0.003798145404262353+0.009953606082753418im
print("relative error = ", (result - ref) / abs(ref), "\n")

## timing
@btime levintegrate($si, $f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
