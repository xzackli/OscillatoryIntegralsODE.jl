
# kernel and functions for the 𝐡𝐚𝐫𝐦𝐨𝐧𝐢𝐜 function exp(i ω x)
struct HarmonicIntegral{T, Tpar} <: OscillatoryIntegral{T, Tpar, 1}
    ω::Tpar
end
HarmonicIntegral{T}(ω::Tpar) where {T, Tpar} = HarmonicIntegral{T,Tpar}(ω)
osc_kernel(oi::HarmonicIntegral{T, Tpar}, x) where {T, Tpar} = 1im * oi.ω
osc_funcs(oi::HarmonicIntegral{T, Tpar}, x) where {T, Tpar} =  # the vector v(x)
    exp(1im * oi.ω * x)


# M=1 case, where the Levin ODE is one-dimensional, like e^iωg(x)
function levintegrate(oi::OscillatoryIntegral{T, Tpar, 1}, f::F, a, b, alg=Vern9();
                      ode_kwargs...) where {T, Tpar, F}
    levin_static(u, p_oi, t) = f(t) - osc_kernel(p_oi, t) * u  # Levin94 eq 2.6
    va = zero(T)  # arbitrary, Levin prescribes no BC
    prob = ODEProblem{false}(levin_static, va, (a, b), oi)  # set up ODE from a to b
    sol = solve(prob, alg, save_start=false, save_everystep=false; ode_kwargs...)
    vb = last(sol)
    return osc_funcs(oi, b) * vb # - osc_funcs(oi, a) * va  # Levin94 eq 2.9, va is 0
end
