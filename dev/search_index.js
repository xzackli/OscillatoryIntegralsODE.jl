var documenterSearchIndex = {"docs":
[{"location":"api_index/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api_index/","page":"API","title":"API","text":"","category":"page"},{"location":"api_index/","page":"API","title":"API","text":"Modules = [OscillatoryIntegralsODE]","category":"page"},{"location":"api_index/#OscillatoryIntegralsODE.levintegrate","page":"API","title":"OscillatoryIntegralsODE.levintegrate","text":"levintegrate(oi, f, a, b, alg=Vern9(); ode_kwargs...)\n\nThis computes int_a^b f(x) w(x) dxfor0 < a < b, wheref(x)is a smooth, univariate function, andw(x)` is a highly-oscillatory function.\n\nArguments:\n\noi::OscillatoryIntegral: parameters of the highly-oscillatory integral\nf: the function to be integrated\na: lower integration bound, must be greater than zero\nb: upper integration bound, must be greater than a\nalg: algorithm used in the ODE solve, use probably want something high order\n\nKeywords\n\nAny keyword arguments are passed to the solve call of the Levin ODE\n\nReturns:\n\nthe result of the integral\n\nExamples\n\njulia> f(x) = exp(-x^2/16)\n       bi = BesselJIntegral{Float64}(100., 100.)  # nu, r\n       levintegrate(bi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)\n\n0.006311599451652101\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = OscillatoryIntegralsODE","category":"page"},{"location":"#OscillatoryIntegralsODE","page":"Home","title":"OscillatoryIntegralsODE","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package computes highly oscillatory integrals by combining a Levin method (Levin 1994) with an ODE solver. It currently supports Bessel, spherical bessel, and harmonic integrals.","category":"page"}]
}
