
# Usage

To calculate a univariate integral with a Bessel function
```math
I = \int_a^b f(x) J_{\nu}(r x) \, dx,
```
one needs to initialize a [`BesselJIntegral`](@ref) with the output type `T`, order ``\nu`` and frequency ``r``. Consider an example,
```math
I = \int_a^b e^{-x^2/16} J_{100}(200 x) \, dx,
```
```julia
using OscillatoryIntegralsODE
# we integrate the Bessel with this f
f(x) = exp(-x^2/16)
# set up problem with output type Float64, Bessel order nu=100, frequency r=200
bi = BesselJIntegral{Float64}(100., 200.)  # nu, r
# now integrate over (a, b) = (1, 5)
levintegrate(bi, f, 1.0, 5.0; abstol=1e-6, reltol=1e-6)
```
```
0.00030145814639532284
```

!!! note
    Levin methods for Bessel functions can only integrate ``0 < a < b``, because ``A(x)`` diverges as it approaches zero. If you need to get down to 0, consider using Gauss quadrature from 0 to the first peak in the Bessel function. Probably the fastest way to do this is to use one of the many asymptotic results for the smallest positive zero (Watson 1922) with ``x_{\nu} = \nu + 1.855757 \nu^{1/3} + O(\nu^{-1/3})``.

We could have used a different ODE integrator. The default [`levintegrate`](@ref) uses `Vern9` from [OrdinaryDiffEq.jl](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Recommended-Methods).

```julia
using OrdinaryDiffEq
levintegrate(bi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
```
Note that the keyword arguments for [`levintegrate`](@ref) are passed directly to the ODE `solve` call.

Similarly for Spherical Bessel functions ``j_{\nu}(rx)``, here's an example integral,
```math
I = \int_1^5 e^{-x^2/16} j_{100}(100 x)\, dx
```

```julia
nu = 100.
r = 100.
f(x) = exp(-x^2/16)
bi = SphericalBesselJIntegral{Float64}(nu, r)
levintegrate(bi, f, 1.0, 5.0; abstol=1e-6, reltol=1e-6)
```
```
0.0008322179291456167
```

We also include an integral type for the harmonic transform ``e^{i \omega x}``,
```julia
ω = 100.
f(x) = exp(-x^2/16)
hi = HarmonicIntegral{Complex{Float64}}(ω)
result = levintegrate(hi, f, 1.0, 5.0; abstol=1e-6, reltol=1e-6)
```
```
0.003797996803011478 + 0.00995319484859839im
```
