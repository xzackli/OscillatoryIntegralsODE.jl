module OscillatoryIntegralsODE

using StaticArrays
using OrdinaryDiffEq
import SpecialFunctions: besselj, sphericalbesselj
import LinearAlgebra: dot

# T is the output type
# Tpar is the type of the parameters (not necessarily the same as output)
# M is the dimension of the Levin ODE, changes for each function
abstract type OscillatoryIntegral{T, Tpar, M} end

"""
    levintegrate(oi, f, a, b, alg=Vern9(); ode_kwargs...)

This computes ``\\int_a^b f(x) w(x) dx`` for ``0 < a < b``, where ``f(x)`` is
a smooth, univariate function, and ``w(x)`` is a highly-oscillatory function.

# Arguments:
- `oi::OscillatoryIntegral`: parameters of the highly-oscillatory integral
- `f`: the function to be integrated
- `a`: lower integration bound, must be greater than zero
- `b`: upper integration bound, must be greater than `a`
- `alg`: algorithm used in the ODE solve, you probably want something high order

# Keywords
- Any keyword arguments are passed to the `solve` call of the Levin ODE

# Returns:
- the result of the integral

# Examples
```julia-repl
julia> f(x) = exp(-x^2/16)
       bi = BesselJIntegral{Float64}(100., 100.)  # nu, r
       levintegrate(bi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)

0.006311599451652101
```
"""
function levintegrate end

# M=2 case, where the Levin ODE is two-dimensional, like Bessel
function levintegrate(oi::OscillatoryIntegral{T, Tpar, 2}, f::F, a, b, alg=Vern9();
                      ode_kwargs...) where {T, Tpar, F}
    levin_static(u, p_oi, t) =
        SA[zero(T); f(t)] .- transpose(osc_kernel(p_oi, t)) * u  # Levin94 eq 2.6
    va = SA[zero(T); zero(T)]  # arbitrary, Levin prescribes no BC
    prob = ODEProblem{false}(levin_static, va, (a, b), oi)  # set up ODE from a to b
    sol = solve(prob, alg, save_start=false, save_everystep=false; ode_kwargs...)
    vb = last(sol)
    return dot(osc_funcs(oi, b), vb) # - dot(osc_funcs(oi, a), va)  # Levin94 eq 2.9, we chose va = 0
end


# OSCILLATORY FUNCTION DEFINITIONS ==================================

# kernel and functions for the 𝐁𝐞𝐬𝐬𝐞𝐥 function J_ν(rx)
struct BesselJIntegral{T, Tpar} <: OscillatoryIntegral{T, Tpar, 2}
    ν::Tpar
    r::Tpar
end
BesselJIntegral{T}(ν::Tpar, r::Tpar) where {T, Tpar} = BesselJIntegral{T,Tpar}(ν, r)
osc_kernel(oi::BesselJIntegral{T, Tpar}, x) where {T, Tpar} =
    SA[(oi.ν-one(Tpar))/x  (-oi.r);
       (oi.r)            -oi.ν/x]  # Levin94 eq 3.3 (the matrix A(x))
osc_funcs(oi::BesselJIntegral{T, Tpar}, x) where {T, Tpar} =
    SA[besselj(oi.ν-one(Tpar), oi.r * x);  besselj(oi.ν, oi.r * x)]  # Levin94 eq 3.1-3.3, v(x)


"""
    BesselJIntegral{T}(ν, r)

This constructs a Bessel integration problem which solves
```math
I = \\int_a^b f(x) J_{\\nu}(r x)\\, dx
```
and produces output of type `{T}`.
"""
function BesselJIntegral end


# kernel and functions for the 𝐬𝐩𝐡𝐞𝐫𝐢𝐜𝐚𝐥 𝐁𝐞𝐬𝐬𝐞𝐥 function j_ν(rx)
struct SphericalBesselJIntegral{T, Tpar} <: OscillatoryIntegral{T, Tpar, 2}
    ν::Tpar
    r::Tpar
end
SphericalBesselJIntegral{T}(ν::Tpar, r::Tpar) where {T, Tpar} =
    SphericalBesselJIntegral{T,Tpar}(ν, r)
osc_kernel(oi::SphericalBesselJIntegral{T, Tpar}, x) where {T, Tpar} =
    SA[((oi.ν-one(Tpar))/x)  (-oi.r);
        (oi.r)            ((-oi.ν-one(Tpar))/x)]  # the matrix A(x)
osc_funcs(oi::SphericalBesselJIntegral{T, Tpar}, x) where {T, Tpar} =  # the vector v(x)
    SA[sphericalbesselj(oi.ν-one(Tpar), oi.r*x);  sphericalbesselj(oi.ν, oi.r*x)]


"""
    SphericalBesselJIntegral{T}(ν, r)

This constructs a spherical Bessel integration problem which solves
```math
I = \\int_a^b f(x) j_{\\nu}(r x)\\, dx
```
and produces output of type `{T}`.

"""
function SphericalBesselJIntegral end


include("harmonic.jl")
export levintegrate, BesselJIntegral, SphericalBesselJIntegral, HarmonicIntegral

end
