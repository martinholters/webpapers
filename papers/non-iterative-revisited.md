---
title: The non-iterative method revisited
author: Martin Holters
---
## Introduction

In [[1]](#ref:Ducceschi-2022), Ducceschi and Bilbao have introduced a family of
non-iterative discretization schemes for ODEs. However, in their presentation,
they make several assumption regarding the structure of the ODE, limiting the
gerenality of these schemes. Here, we show that for the second-order accurate
variant in particular, the method is applicable to general (potentially
vector-valued) ODEs.

In the following, `$x(t)$` denotes the state of a system goverened by an ODE and
`$\dot{x}(t)$` it derivative. When obvious from the context, we drop the
argument `$t$`. We seek an approximation `$\hat{x}(n) \approx x(nT)$` where
`$T$` is the sampling interval. For given signals like an input `$u(t)$`, we
likewise denote samples as e.g. `$\hat{u}(n)=u(nT)$`.

## Recap of the non-iterative method
We start with a recap of the main take-aways of [[1]](#ref:Ducceschi-2022).
Consider a scalar ODE
```equation label=ode
\dot{x} + f(x) = u
```
where `$f(0)=0$`, so that we can rewrite as
```equation
\dot{x} + g(x)\cdot x = u
```
by letting `$g(x)=f(x)/x$` (with continuous extension at `$x=0$`). For the
zero-input case, the non-iterative scheme is then given by
```equation label=non-iterative-implicit
\sigma_P(\hat{x}(n))\cdot\frac{1}{T}(\hat{x}(n+1)-\hat{x}(n))
+ g(\hat{x}(n))\cdot\frac{1}{2}(\hat{x}(n+1)+\hat{x}(n)) = 0
```
where
```equation
\sigma_P(\hat{x}(n)) = \sum_{p=0}^P T^p \zeta_p(\hat{x}(n))
```
and appropriate choice of `$\zeta_p$` leads to a method of accuracy order
`$P+1$`. In particular `$\zeta_0(\hat{x}(n))=1$` and
```equation
\zeta_1(\hat{x}(n))=\frac{1}{2}\left(f'(\hat{x}(n))-g(\hat{x}(n))\right)
```
where `$f'=\frac{df}{dx}$`.

To give the non-iterative scheme, {% eqref "non-iterative-implicit" %} is solved
as
```equation
\hat{x}(n+1)
= \frac{\sigma_P(\hat{x}(n)) - \frac{T}{2}g(\hat{x}(n))}{\sigma_P(\hat{x}(n))
+ \frac{T}{2}g(\hat{x}(n))}\cdot\hat{x}(n)
= \frac{1-\kappa(n)}{1+\kappa(n)}\cdot\hat{x}(n)
```
with
```equation
\kappa(n) = \frac{Tg(\hat{x}(n))}{2\sigma_P(\hat{x}(n))}.
```
For non-zero input, {% eqref "non-iterative-implicit" %} is extended to
```equation
\begin{split}
\sigma_P(\hat{x}(n))\cdot\frac{1}{T}(\hat{x}(n+1)-\hat{x}(n))
+ g(\hat{x}(n))\cdot\frac{1}{2}(\hat{x}(n+1)+\hat{x}(n)) \quad \\
= \frac{1}{2}(\hat{u}(n+1)+\hat{u}(n))
\end{split}
```
which likewise can be solved to
```equation
\begin{split}
\hat{x}(n+1)
&=
\frac{\sigma_P(\hat{x}(n))-\frac{T}{2}g(\hat{x}(n))}
     {\sigma_P(\hat{x}(n))+ \frac{T}{2}g(\hat{x}(n))}
\hat{x}(n)\\
&\quad + \frac{\frac{T}{2}}{\sigma_P(\hat{x}(n))+ \frac{T}{2}g(\hat{x}(n))}
  (\hat{u}(n+1)+\hat{u}(n)).
\end{split}
```

For the vector case, consider an ODE like {% eqref "ode" %} with
```equation
f(x) = Bx + Fq(\eta(x))
```
where `$\eta(x) = F^Tx + c$` where `$c$` is a time-dependent input function and
`$q(\eta) = \begin{pmatrix} q_1(\eta_1) & \cdots & q_N(\eta_N) \end{pmatrix}^T$`
is an element-wise non-linear function. Then in {% eqref
"non-iterative-implicit" %}, we use
```equation
g(x) = B+FD(x)F^T
```
where `$D(x) = \operatorname{diag}(q_i(\eta_i)/\eta_i)$` with
`$\eta(\hat{x}(n)) = F^T\hat{x}(n) + \frac{1}{2}(\hat{c}(n+1)+\hat{c}(n))$` and
`$f'$` becomes the Jacobian.

## Revisiting the second-order case

We now only consider an accuracy order of `$2$`, we thus need `$P=1$`, i.e.
```equation
\sigma_1(\hat{x}(n)) = 1 + \frac{T}{2}\left(f'(\hat{x}(n))-g(\hat{x}(n))\right).
```
We again consider an ODE in the form of {% eqref "ode" %}, but vector-valued and
assume `$g(x)x=f(x)$`, where matrix `$g(x)$` is not uniquely defined. Focusing
on the zero-input case first, we rewrite {% eqref "non-iterative-implicit" %} as
```equation 
\left(\sigma_1(\hat{x}(n)) + \frac{T}{2}g(\hat{x}(n))\right)\hat{x}(n+1)
=
\left(\sigma_1(\hat{x}(n))-\frac{T}{2}g(\hat{x}(n))\right)\hat{x}(n)
```
and substitute `$\sigma_1$` to obtain
```equation 
\left(I + \frac{T}{2}f'(\hat{x}(n))\right)\hat{x}(n+1)
=
\left(I + \frac{T}{2}f'(\hat{x}(n))-Tg(\hat{x}(n))\right)\hat{x}(n).
```
We can now exploit `$g(x)x=f(x)$` to rewrite to
```equation 
\left(I + \frac{T}{2}f'(\hat{x}(n))\right)\hat{x}(n+1)
=
\left(I + \frac{T}{2}f'(\hat{x}(n))\right)\hat{x}(n)-Tf(\hat{x}(n)).
```
and solve as
```equation label=non-iterative-order-2-explicit
\hat{x}(n+1)
=
\hat{x}(n)-\left(I + \frac{T}{2}f'(\hat{x}(n))\right)^{-1}Tf(\hat{x}(n)).
```
We have thus obtained a non-iterative scheme without explicit appearance of
`$g$`, which even lets us treat cases where `$f(0)\ne0$`. It remains to be
verified that the scheme is second-order accurate even then.

For sufficiently small `$T$`, we may utilize the geometric series
```equation
\sum_{i=0}^\infty \left(-\tfrac{T}{2}f'(\hat{x}(n))\right)^i
= \left(I + \frac{T}{2}f'(\hat{x}(n))\right)^{-1}
```
so that
```equation 
\begin{split}
\hat{x}(n+1)
&=
\hat{x}(n)
- \left(
    \sum_{i=0}^\infty\left(-\tfrac{T}{2}f'(\hat{x}(n))\right)^i
  \right)Tf(\hat{x}(n)) \\
&= \hat{x}(n)-Tf(\hat{x}(n))
   + \frac{T^2}{2}f'(\hat{x}(n))f(\hat{x}(n)) + O(T^3).
\end{split}
```
Comparing with the Taylor series expansion
```equation
\begin{split}
x(t+T) &= x(t) + T\dot{x}(t) + \frac{T^2}{2}\ddot{x}(t) + O(T^3) \\
&= x(t) - Tf(x(t)) + \frac{T^2}{2}f'(x(t))f(x(t)) + O(T^3)
\end{split}
```
it is obvious that the method is second-order accurate.

### Including a driving term
Now consider a vector-valued ODE of the form
```equation
\dot{x} + f(x,c) = 0
```
where `$f$` is arbitrary and includes a driving term `$c$`. (For simplicity, we
only consider a zero right-hand side, as a non-zero right-hand side could easily
be subsumed in `$f$`.) In the following, let `$J_x$` and `$J_c$` denote the
Jacobians of `$f$` with respect to `$x$` and `$c$`, repectively. We adapt
{% eqref "non-iterative-order-2-explicit" %} as
```equation 
\hat{x}(n+1)
=
\hat{x}(n)-\left(I + \frac{T}{2}J_x^n\right)^{-1}
Tf^n
```
where `$f^n=f\big(\hat{x}(n), \tfrac{1}{2}(\hat{c}(n+1)+\hat{c}(n))\big)$` and
likewise
`$J_x^n=J_x\big(\hat{x}(n),\tfrac{1}{2}(\hat{c}(n+1)+\hat{c}(n))\big)$`. By
following the same reasoning as above, we find
```equation label=non-iterative-order-2-taylor
\hat{x}(n+1)
=
\hat{x}(n)-Tf^n + \frac{T^2}{2}J_x^nf^n + O(T^3).
```
The Taylor series of `$x(t+T)$` becomes a bit more complicated due to the second
(time-dependent) argument to `$f$`. In particular. we get
```equation
\begin{split}
x(t+T) &= x(t) + T\dot{x}(t) + \frac{T^2}{2}\ddot{x}(t) + O(T^3) \\
&= x(t) - Tf(x(t),c(t)) \\
& \quad + \frac{T^2}{2}\big(
    J_x(x(t),c(t))f(x(t),c(t)) - J_c(x(t),c(t))\dot{c}(t)
  \big)
  + O(T^3).
\end{split}
```
We now utilize the Taylor expansion
```equation
f(x(t),\tfrac12(c(t+T)+c(t)))
= f(x(t),c(t)) + \frac{T}{2}J_c(x(t),c(t))\dot{c}(t) + O(T^2)
```
to replace
```equation
f(x(t),c(t))
=
f(x(t),\tfrac12(c(t+T)+c(t)))
- \frac{T}{2}J_c(x(t),c(t))\dot{c}(t) + O(T^2)
```
yielding
```equation
\begin{split}
x(t+T)
&= x(t) - Tf(x(t),\tfrac12(c(t+T)+c(t))) \\
& \quad + \frac{T^2}{2}J_x(x(t),c(t))f(x(t),\tfrac12(c(t+T)+c(t)))+ O(T^3).
\end{split}
```
Note that the terms involving `$J_c$` have cancelled and we now have `$f$`
appearing with the desired arguments. We note in passing that if we had used
`$J_x^n=J_x\big(\hat{x}(n), \hat{c}(n)\big)$`, this would already verify
second-order accurateness. However, we have to do one more step and employ the
Taylor series
```equation
J_x(x(t),\tfrac12(c(t+T)+c(t))) = J_x(x(t),c(t)) + O(T)
```
to further rewrite to
```equation
\begin{split}
x(t+T)
&= x(t) - Tf(x(t),\tfrac12(c(t+T)+c(t))) \\
& \quad + \frac{T^2}{2}J_x(x(t),\tfrac12(c(t+T)+c(t)))
          f(x(t),\tfrac12(c(t+T)+c(t))) \\
& \quad + O(T^3).
\end{split}
```
Comparison with {% eqref "non-iterative-order-2-taylor" %} now shows that the
method is indeed second-order accurate.

## Bibliography

<div id="ref:Ducceschi-2022">[1] M. Ducceschi and S. Bilbao,
Non-iterative simulation methods for virtual analog modelling.
<i>IEEE/ACM Transactions on Audio, Speech, and Language Processing</i>,
vol. 50, pp. 3189–3198, 2022.</div>
