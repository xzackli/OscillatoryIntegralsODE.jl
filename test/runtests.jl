using OscillatoryIntegralsODE
using OrdinaryDiffEq
using Test

@testset "bessel" begin
    nu = 100.
    r = 100.
    f(x) = exp(-x^2/16)
    bi = BesselJIntegral{Float64}(nu, r)
    result = levintegrate(bi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
    ref = 0.006311630277650583
    @test abs( (result - ref) / ref ) < 1e-4
end

@testset "spherical bessel" begin
    nu = 100.
    r = 100.
    f(x) = exp(-x^2/16)
    bi = SphericalBesselJIntegral{Float64}(nu, r)
    result = levintegrate(bi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
    ref = 0.0008322179291456167
    @test abs( (result - ref) / ref ) < 1e-4
end

@testset "harmonic" begin
    ω = 100.
    f(x) = exp(-x^2/16)
    hi = HarmonicIntegral{Complex{Float64}}(ω)
    result = levintegrate(hi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
    ref = 0.003798145404262353+0.009953606082753418im
    @test abs(result - ref) / abs(ref) < 1e-4
end
