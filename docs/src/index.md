```@meta
CurrentModule = OscillatoryIntegralsODE
```

# OscillatoryIntegralsODE

This package performs highly oscillatory quadrature by combining a Levin method ([Levin 1994](https://www.sciencedirect.com/science/article/pii/0377042794001189)) with an ODE solver. It currently supports Bessel, spherical Bessel, and harmonic integrals. 

## Background

Consider the commonly-encountered class of integrals

```math
I = \int_a^b g(x) S(rx) \, dx
```

where ``g(x)`` is a well-behaved function, but ``S`` is rapidly-oscillating for large values of ``r``. In cosmology, the spherical Bessel function ``j_{\nu}(rx)`` is a common example of a highly oscillatory ``S``. Traditional Gaussian quadrature has poor convergence properties for such integrals, as the number of required nodes scales with frequency ``r``. There is a rich literature on other quadrature approaches (for a review, see [Huybrechs and Olver 2012](https://people.cs.kuleuven.be/~daan.huybrechs/research/HOQ.pdf)). This package uses the Levin method in combination with an ODE solver to efficiently evaluate integrals with ``S`` functions of Bessel ``J_{\nu}(rx)``, spherical Bessel ``j_{\nu}(rx)``, and harmonic ``e^{i \omega x}`` for finite ``0 < a < b``.

The main idea of the Levin method ([Levin 1994](https://www.sciencedirect.com/science/article/pii/0377042794001189)) is to transform the integral of a rapidly-oscillating function into an equivalent non-rapidly oscillatory system of ordinary differential equations. We now follow Levin, and let ``\mathbf{f}(x) = ( f_1(x), \ldots, f_m(x) )^T`` be an ``m``-vector of non-rapidly oscillatory functions, and ``\mathbf{w}(x) = ( w_1(x), \ldots, w_m(x) )^T`` be an ``m``-vector of linearly independent rapidly oscillatory functions. Now suppose your oscillatory functions ``w_i`` have some special properties, such that the derivative vector ``\mathbf{w}'(x)`` is always some linear transformation of ``mathbf{w}``. Then we could write 

```math
\mathbf{w}'(x) = \mathbf{A}(x) \mathbf{w}(x).
```
We are lucky, as many common rapidly oscillatory functions, such as those of the Bessel, spherical Bessel, harmonic, sine, and cosine transforms, satisfy such a relation. For example, the derivative of the Bessel function can be expressed as a linear combination of two Bessel functions, and they are also subject to a convenient recursion relation among different orders. One only needs ``m=2`` to be able to write a ``\mathbf{w}(x)`` containing a desired Bessel function (see later sections).

Now suppose there exists an ``m``-vector of functions ``\mathbf{p}(x)`` such that 
```math
\frac{d}{dx} \left( \mathbf{p} \cdot \mathbf{w} \right) \approx \mathbf{f} \cdot \mathbf{w}.
```
If we could obtain such a ``\mathbf{p}``, then by the fundamental theorem of calculus, our integral becomes
```math
I \approx \int_a^b \frac{d}{dx} \left( \mathbf{p}   \cdot \mathbf{w} \right) = \mathbf{p}(b) \cdot \mathbf{b} - \mathbf{p}(a) \cdot \mathbf{a}.
```
We have thus reduced our integral to finding an approximation to functions ``\mathbf{p}(x)``. Applying the product rule, we have (recall that we can write ``\mathbf{w}^\prime = \mathbf{A} \mathbf{w}``)
```math
(\mathbf{p} \cdot \mathbf{w})^\prime = \mathbf{p}^\prime \cdot \mathbf{w} + \mathbf{p} \cdot \mathbf{w}^\prime = \mathbf{p}^\prime \cdot \mathbf{w} + \mathbf{p} \cdot \mathbf{A} \mathbf{w} = (\mathbf{p}^\prime + \mathbf{A}^T \mathbf{p}) \cdot \mathbf{w} \approx \mathbf{f} \cdot \mathbf{w}.
```
Thus ``\mathbf{p}`` should be an approximate solution to 
```math
\mathbf{p}^\prime + \mathbf{A}^T \mathbf{p} = \mathbf{f}.
```
We solve this ODE using an ODE solver. As there are no boundary conditions here, we generally choose the initial condition ``p(a) = \mathbf{0}``. We then only need to obtain the ODE solution at ``b``, and we have our integral ``I \approx \mathbf{p}(b) \cdot \mathbf{w}(b)``. If the ODE solver is a Runge Kutta method, then there's a good chance it corresponds to some kind of Taylor-based collocation method, which is typically what the literature discuses.

### Bessel Functions
Bessel functions of the first kind ``J_{\nu}`` satisfy 
```math
J^{\prime}_{\nu - 1}(x) = \frac{\nu - 1}{x} J_{\nu - 1}(x) - J_{\nu}(x), \quad \quad
J^{\prime}_{\nu}(x) =  J_{\nu - 1}(x) - \frac{\nu}{x} J_{\nu}(x)
```
That means we can have ``\mathbf{w}(x) = (J_{\nu - 1}(x), J_{\nu}(x))`` and ``\mathbf{w}^\prime = \mathbf{A} \mathbf{w}`` if 
```math
\mathbf{A}(x) = \begin{pmatrix}
(\nu - 1) / x & -r \\ 
r &  - \nu / x
\end{pmatrix}
```
We define ``\mathbf{f}(x) = (0, f(x))``. We now have our two variables ``\mathbf{p}^\prime + \mathbf{A}^T \mathbf{p} = \mathbf{f}``. Once we solve that ODE, we can get our result ``I \approx \mathbf{p}(b) \cdot \mathbf{w}(b)``.

### Spherical Bessel Functions
Spherical Bessel functions of the first kind ``j_{\nu}`` satisfy a similar differential relation, such that for ``\mathbf{w}(x) = (j_{\nu - 1}(x), j_{\nu}(x))``, we have
```math
\mathbf{A}(x) = \begin{pmatrix}
(\nu - 1) / x & -r \\ 
r &  - (\nu+1) / x
\end{pmatrix}
```

### Harmonic Integrals
The function ``e^{i \omega x}`` requires only ``m=1``. In this case we have a scalar ``A``,
```math
A(x) = i \omega
```
for ``w(x) = e^{i \omega x}``, and pretty trivially ``w^\prime(x) = A(x) w(x)``.

## Usage

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
# set up problem with output type Float64, Bessel order nu=100, frequency r=100
bi = BesselJIntegral{Float64}(100., 200.)  # nu, r
# now integrate over (a, b) = (1, 5)
levintegrate(bi, f, 1.0, 5.0; abstol=1e-6, reltol=1e-6)
```
```
0.00002612881708428357
```

!!! note
    Levin methods for Bessel functions can only integrate ``0 < a < b``, because ``A(x)`` diverges as it approaches zero. If you need to get down to 0, consider using Gauss quadrature from 0 to the first peak in the Bessel function. Probably the fastest way to do this is to use one of the many asymptotic results for the smallest positive zero, i.e. Watson 1922 ``x_{\nu} = \nu + 1.855757 \nu^{1/3} + O(\nu^{-1/3})``.

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
levintegrate(bi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
```
```
0.0008322179291456167
```

We also include an integral type for the harmonic transform ``e^{i \omega x}``,
```julia
ω = 100.
f(x) = exp(-x^2/16)
hi = HarmonicIntegral{Complex{Float64}}(ω)
result = levintegrate(hi, f, 1.0, 5.0, Vern9(); abstol=1e-6, reltol=1e-6)
```
```
0.003797996803011478 + 0.00995319484859839im
```
