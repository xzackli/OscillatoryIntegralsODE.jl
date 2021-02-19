# OscillatoryIntegralsODE

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xzackli.github.io/OscillatoryIntegralsODE.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xzackli.github.io/OscillatoryIntegralsODE.jl/dev)
[![Build Status](https://github.com/xzackli/OscillatoryIntegralsODE.jl/workflows/CI/badge.svg)](https://github.com/xzackli/OscillatoryIntegralsODE.jl/actions)
[![Coverage](https://codecov.io/gh/xzackli/OscillatoryIntegralsODE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xzackli/OscillatoryIntegralsODE.jl)

This package computes highly oscillatory integrals by combining a Levin method ([Levin 1994](https://www.sciencedirect.com/science/article/pii/0377042794001189)) with an ODE solver. It currently supports Bessel, spherical bessel, and harmonic integrals.

```julia
using OscillatoryIntegralsODE
f(x) = exp(-x^2/16)
bi = BesselJIntegral{Float64}(100., 100.)  # nu, r
levintegrate(bi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
```
```
0.006311599451652101
```