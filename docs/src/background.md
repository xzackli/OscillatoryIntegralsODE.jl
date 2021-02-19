
# Background

Consider the commonly-encountered class of integrals

```math
I = \int_a^b g(x) S(rx) \, dx
```

where ``g(x)`` is a well-behaved function, but ``S`` is rapidly-oscillating for large values of ``r``. In cosmology, the spherical Bessel function ``j_{\nu}(rx)`` is a common example of a highly oscillatory ``S``. Traditional Gaussian quadrature has poor convergence properties for such integrals, as the number of required nodes scales with frequency ``r``. There is a rich literature on other quadrature approaches (for a review, see [Huybrechs and Olver 2012](https://people.cs.kuleuven.be/~daan.huybrechs/research/HOQ.pdf)). This package uses the Levin method in combination with an ODE solver to efficiently evaluate integrals with where ``S`` is a Bessel ``J_{\nu}(rx)``, spherical Bessel ``j_{\nu}(rx)``, or harmonic transform ``e^{i \omega x}`` for finite ``0 < a < b``.

The main idea of the Levin method ([Levin 1994](https://www.sciencedirect.com/science/article/pii/0377042794001189)) is to transform the integral of a rapidly-oscillating function into an equivalent non-rapidly oscillatory system of ordinary differential equations. We now follow Levin, and let ``\mathbf{f}(x) = ( f_1(x), \ldots, f_m(x) )^T`` be an ``m``-vector of non-rapidly oscillatory functions, and ``\mathbf{w}(x) = ( w_1(x), \ldots, w_m(x) )^T`` be an ``m``-vector of linearly independent rapidly oscillatory functions. Now suppose your oscillatory functions ``w_i`` have some special properties, such that the derivative vector ``\mathbf{w}'(x)`` is always some linear transformation of ``\mathbf{w}``. Then we could write 

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
We solve this ODE using an ODE solver. As there are no boundary conditions here, we generally choose the initial condition ``p(a) = \mathbf{0}``. We then only need to obtain the ODE solution at ``b``, and we have our integral ``I \approx \mathbf{p}(b) \cdot \mathbf{w}(b)``. If the ODE solver is a Runge-Kutta method, then there's a good chance it corresponds to some kind of Taylor-based collocation method, which is typically what the literature discuses.

### Bessel Functions
Bessel functions of the first kind ``J_{\nu}`` satisfy 
```math
J^{\prime}_{\nu - 1}(x) = \frac{\nu - 1}{x} J_{\nu - 1}(x) - J_{\nu}(x), \quad \quad
J^{\prime}_{\nu}(x) =  J_{\nu - 1}(x) - \frac{\nu}{x} J_{\nu}(x)
```
That means we can have ``\mathbf{w}(x) = (J_{\nu - 1}(r x), J_{\nu}(r x))`` and ``\mathbf{w}^\prime = \mathbf{A} \mathbf{w}`` if 
```math
\mathbf{A}(x) = \begin{pmatrix}
(\nu - 1) / x & -r \\ 
r &  - \nu / x
\end{pmatrix}
```
We define ``\mathbf{f}(x) = (0, f(x))``. We now have our two variables ``\mathbf{p}^\prime + \mathbf{A}^T \mathbf{p} = \mathbf{f}``. Once we solve that ODE, we can get our result ``I \approx \mathbf{p}(b) \cdot \mathbf{w}(b)``.

### Spherical Bessel Functions
Spherical Bessel functions of the first kind ``j_{\nu}`` satisfy a similar differential relation, such that for ``\mathbf{w}(x) = (j_{\nu - 1}(r x), j_{\nu}(r x))``, we have
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
