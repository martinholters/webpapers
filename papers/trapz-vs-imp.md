---
title: Trapezoidal rule vs. implicit midpoint method
author: Martin Holters
---

## Introduction

Consider an autonomous ordinary differetial equation (ODE)
```equation
\dot{x} = f(x)
```
where `$x$` is the state vector. For general `$f$`, sovling analytically is
usually impossible and we have to resort to numerical methods that seek to
compute a sequence `$x^n\approx x(nT)$`, where `$T$` is the time-step, here
assumed constant for simplicity. Furthermore, we only consider the initial value
problem where `$x(0)=x_0$` is given and use `$x^0 = x_0$`.

There are many methods to achieve this goal of which two will be considered
here: the trapezoidal rule and the implicit midpoint method. In particular, we
shall review their similarity.

## Reminder of the implicit midpoint rule

Just as a brief reminder, for an ODE `$\dot{x} = f(x)$`, the implcit midpoint
method is given by
```equation label=imp
  x^{n+1} = x^n + Tf\big(\tfrac{1}{2}(x^{n+1} + x^n)\big).
```

## The trapezoidal rule revisited

The trapezoidal rule applied it to an ODE `$\dot{x} = f(x)$` is usually stated as
```equation label=tr
  x^{n+1} = x^n + \tfrac{T}{2}\left(f(x^{n+1}) + f(x^n)\right).
```
But for some applications, another form is more suitable which will be derived
next. We start by rewriting {% eqref "tr" %} as
```equation label=xhalf
  x^{n+1} - \tfrac{T}{2}f(x^{n+1}) = x^n + \tfrac{T}{2}f(x^n) =: x^{n+\frac12}.
```
In passing, we note that this gives rise to an interesting interpretation of
the trapezoidal rule: A forward-Euler step of length `$\frac{T}{2}$` to reach
`$x^{n+\frac12}$` followed by a backward-Euler step also of length `$\frac{T}{2}$`
to reach `$x^{n+1}$`. But here, we use it to rewrite the trapezoidal rule in
`$x^{n+\frac12}$` and `$x^{n-\frac12}$`, where the latter is, from the left-hand side of
{% eqref "xhalf" %}, given by
```equation
  x^{n-\frac12} = x^{n} - \tfrac{T}{2}f(x^{n}).
```
Combining this with the right-hand side definition of {% eqref "xhalf" %} for
`$x^{n+\frac12}$`, we therefore get
```equation
\begin{split}
  x^{n+\frac12} + x^{n-\frac12} &= 2x^n \\
  x^{n+\frac12} - x^{n-\frac12} &= Tf(x^{n})
\end{split}
```
or equivalently
```equation label=trapzstate
  x^n = \tfrac12 \left(x^{n+\frac12} + x^{n-\frac12}\right)
```
and
```equation label=trapzf
  f(x^{n}) = \tfrac{1}{T}\left(x^{n+\frac12} - x^{n-\frac12}\right).
```
Inserting {% eqref "trapzstate" %} in {% eqref "trapzf" %} we get the
alternative trapezoidal rule formulation
```equation
  \tfrac{1}{T}\left(x^{n+\frac12} - x^{n-\frac12}\right)
  = f\left(\tfrac12 \left(x^{n+\frac12} + x^{n-\frac12}\right)\right)
```
or
```equation label=trapzalt
  x^{n+\frac12}
  = x^{n-\frac12} + Tf\left(\tfrac12 \left(x^{n+\frac12} + x^{n-\frac12}\right)\right).
```
This is very similar to the implicit midpoint rule, but with a half-sample time
shift. Also note that from the sequence obtained from {% eqref "trapzalt" %},
the actual states have to be recovered with {% eqref "trapzstate" %}.

## The implicit midpoint method revisited

Given that we could rewrite the trapezoidal rule such that it looked like the
implicit midpoint rule time-shifted by `$\frac{1}{2}$`, it should be possbile to
also go the other way. Indeed, letting
```equation label=impstate
x^{n+\tfrac12} = \frac{1}{2}(x^{n+1}+x^{n})
```
we can rewrite {% eqref "imp" %} as
```equation label=impupdate1
  x^{n+1} = x^n + Tf\big(x^{n+\tfrac12}\big)
```
i.e.
```equation label=impforw
  Tf\big(x^{n+\tfrac12}\big) = x^{n+1} - x^n = 2(x^{n+\tfrac12} - x^n)
```
utilizing `$x^{n+1}=2x^{n+\tfrac12}-x^{n}$` obtained from {% eqref "impstate"
%}. Time-shifting {% eqref "impupdate1" %} and utilizing
`$x^{n-1}=2x^{n-\tfrac12}-x^{n}$` from a time-shifted {% eqref "impstate" %}, we
similarly get
```equation label=impbackw
  Tf\big(x^{n-\tfrac12}\big) = (x^{n} - x^{n-1}) = 2(x^{n} - x^{n-\tfrac12}).
```
Summing {% eqref "impforw" %} and {% eqref "impbackw" %} yields
```equation
  Tf\big(x^{n+\tfrac12}\big) + Tf\big(x^{n-\tfrac12}\big)
  = 2(x^{n+\tfrac12} - x^{n-\tfrac12})
```
which can be rewritten to
```equation label=impalt
  x^{n+\tfrac12}
  = x^{n-\tfrac12}
    + \frac{T}{2}\left(f\big(x^{n+\tfrac12}\big) + f\big(x^{n-\tfrac12}\big)\right),
```
indeed resembling a time-shifted version of the trapezoidal rule. To recover
`$x^n$`, we can solve {% eqref "impforw" %} or {% eqref "impbackw" %}
to
```equation label=imp_rec1
  x^n = x^{n+\tfrac12} - \frac{T}{2}f\big(x^{n+\tfrac12}\big)
```
or
```equation label=imp_rec2
  x^{n} = x^{n-\tfrac12} + \frac{T}{2}f\big(x^{n-\tfrac12}\big),
```
respectively.

## Example

Consider
```equation
f(x) =- x^3.
```
Then the trapezoidal rule in traditional form reads
```equation label=example_trapz
x_\text{TR}^{n+1} = x_\text{TR}^n + \frac{T}{2}\big(f(x_\text{TR}^{n+1})+f(x_\text{TR}^{n})\big)
= x_\text{TR}^n - \frac{T}{2}\big((x_\text{TR}^{n+1})^3+(x_\text{TR}^{n})^3\big).
```
For the alternative formulation, we get
```equation label=example_trapzalt
x_\text{TR}^{n+\frac12} = x_\text{TR}^{n-\frac12} + Tf\big(\tfrac{1}{2}(x_\text{TR}^{n+\frac12}+x_\text{TR}^{n-\frac12})\big)
= x_\text{TR}^{n-\frac12} - T\big(\tfrac{1}{2}(x_\text{TR}^{n+\frac12}+x_\text{TR}^{n-\frac12})\big)^3.
```
But we need to recover the original states with
`$x_\text{TR}^{n}=\frac{1}{2}\big(x_\text{TR}^{n+\frac12}+x_\text{TR}^{n-\frac12}\big)$`
or equivalently
`$x_\text{TR}^{n+1}=\frac{1}{2}\big(x_\text{TR}^{n+\frac32}+x_\text{TR}^{n+\frac12}\big)$`
which expands to
```equation
\begin{split}
x_\text{TR}^{n+1}&=\frac{1}{2}\big(x_\text{TR}^{n+\frac32}+x_\text{TR}^{n+\frac12}\big) \\
&= \frac{1}{2}\cdot\left(x_\text{TR}^{n+\frac12} - T\big(\tfrac{1}{2}(x_\text{TR}^{n+\frac32}+x_\text{TR}^{n+\frac12})\big)^3 + x_\text{TR}^{n+\frac12}\right) \\
&= x_\text{TR}^{n+\frac12} - \frac{T}{2}\big(\tfrac{1}{2}(x_\text{TR}^{n+\frac32}+x_\text{TR}^{n+\frac12})\big)^3 \\
&= x_\text{TR}^{n+\frac12} - \frac{T}{2}\big(x_\text{TR}^{n+1}\big)^3
\end{split}
```
where the same identity has been used in the last step. Now remember that
`$x_\text{TR}^{n+\tfrac12}=x_\text{TR}^n+\frac{T}{2}f(x_\text{TR}^n)=x_\text{TR}^n-\frac{T}{2}(x_\text{TR}^n)^3$`
and we recover
`$x_\text{TR}^{n+1} = x_\text{TR}^n-\frac{T}{2}(x_\text{TR}^n)^3-\frac{T}{2}(x_\text{TR}^{n+1})^3$`,
i.e. the same output as for the traditional variant as expected.

In contrast, the implicit midpoint method yields
```equation label=example_imp
x_\text{IMP}^{n+1} = x_\text{IMP}^n + Tf\big(\tfrac{1}{2}(x_\text{IMP}^{n+1}+x_\text{IMP}^{n})\big)
= x_\text{IMP}^n - T\big(\tfrac{1}{2}(x_\text{IMP}^{n+1}+x_\text{IMP}^{n})\big)^3
```
and for the alternative formulation
```equation label=example_impalt
x_\text{IMP}^{n+\tfrac12} = x_\text{IMP}^{n-\tfrac12} + \frac{T}{2}\big(f(x_\text{IMP}^{n+\tfrac12})+f(x_\text{IMP}^{n-\tfrac12})\big)
= x_\text{IMP}^{n-\tfrac12} - \frac{T}{2}\big((x_\text{IMP}^{n+\tfrac12})^3+(x_\text{IMP}^{n-\tfrac12})^3\big).
```
Recovering `$x_\text{IMP}^{n+1}$` from the latter with {% eqref "imp_rec2" %}, we get
```equation
\begin{split}
x_\text{IMP}^{n+1}
&= x_\text{IMP}^{n+\tfrac12} + \frac{T}{2}f\big(x_\text{IMP}^{n+\tfrac12}\big) \\
&= x_\text{IMP}^{n+\tfrac12} - \frac{T}{2}\big(x_\text{IMP}^{n+\tfrac12}\big)^3 \\
\end{split}
```
substituting {% eqref "impstate" %} the gives
```equation
x_\text{IMP}^{n+1}
= \frac{1}{2}(x_\text{IMP}^{n+1}+x_\text{IMP}^{n}) - \frac{T}{2}\big(\tfrac{1}{2}(x_\text{IMP}^{n+1}+x_\text{IMP}^{n})\big)^3
```
which can be solved to equal {% eqref "example_imp" %}.


The following table shows the results for the approaches with `$T=1$` starting
from `$x_\text{TR}^0=x_\text{IMP}^0=x_0=1$`,
`$x_\text{TR}^{\tfrac12} = x_0+\frac{1}{2}(-x_0)^3=0.5$`, and
`$x_\text{IMP}^{\tfrac12} = \frac{1}{2}(x_0+x_\text{IMP}^1)=0.77092$`:

`$n$` | `$x_\text{TR}^n$` | `$x_\text{TR}^{n+\tfrac12}$` | `$x_\text{IMP}^n$` | `$x_\text{IMP}^{n+\tfrac12}$`
------|-------------------|------------------------------|--------------------|------------------------------
{% for row in example_results
%}{{
  row.n
}} | {{
  row.trapz | round: 5
}} | {{
  row.trapzalt | round: 5
}} | {{
  row.imp | round: 5
}} | {{
  row.impalt | round: 5
}}
{% endfor %}

Here, the `$x_\text{TR}^n$` column is computed using
{% eqref "example_trapz" %}, the `$x_\text{TR}^{n+\tfrac12}$` column is computed
using {% eqref "example_trapzalt" %}, the `$x_\text{IMP}^n$` column is
computed using {% eqref "example_imp" %}, and the `$x_\text{IMP}^{n+\tfrac12}$`
column is computed using {% eqref "example_impalt" %}. It can be verified that
indeed
`$x_\text{TR}^n=\frac{1}{2}(x_\text{TR}^{n+\tfrac12}+x_\text{TR}^{n-\tfrac12})$`
and likewise for the implicit midpoint method, we may verify that
`$x_\text{IMP}^{n+\tfrac12}=\frac{1}{2}(x_\text{IMP}^{n+1}+x_\text{IMP}^{n})$`.
This can also be seen in the following figure: For the trapezoidal rule, the
`$x_\text{TR}^n$` points are located on line segments joining
`$x_\text{TR}^{n\pm\tfrac12}$`, while for the implicit midpoint method, it is
the other way round.
<svg version="1.1" width="310" height="245" class="figure" xmlns="http://www.w3.org/2000/svg">
  <defs>
    <marker id="circle" markerWidth="8" markerHeight="8" refX="4" refY="4">
      <circle cx="4" cy="4" r="3.5" stroke="#000" fill="none" />
    </marker>
    <marker id="disc" markerWidth="8" markerHeight="8" refX="4" refY="4">
      <circle cx="4" cy="4" r="4" stroke="none" fill="#000" />
    </marker>
    <marker id="x" markerWidth="8" markerHeight="8" refX="4" refY="4">
      <line x1="0" y1="0" x2="8" y2="8" stroke="#000" />
      <line x1="0" y1="8" x2="8" y2="0" stroke="#000" />
    </marker>
    <marker id="+" markerWidth="8" markerHeight="8" refX="4" refY="4">
      <line x1="0" y1="4" x2="8" y2="4" stroke="#000" />
      <line x1="4" y1="0" x2="4" y2="8" stroke="#000" />
    </marker>
  </defs>
  <g transform="translate(30, 10)">
    <polyline stroke="none" fill="none" marker-start="url(#x)" marker-mid="url(#x)" marker-end="url(#x)" points="
      {% for row in example_results %}{{
        row.n | times: 50
      }} {{
        row.trapz | times: -200 | plus: 200
      }}
      {% endfor %}" />
    <polyline stroke="#000" fill="none" marker-start="url(#circle)" marker-mid="url(#circle)" marker-end="url(#circle)" points="
      {% for row in example_results %}{{
        row.n | plus: 0.5 | times: 50
      }} {{
        row.trapzalt | times: -200 | plus: 200
      }}
      {% endfor %}" />
    <polyline stroke="#000" fill="none" marker-start="url(#disc)" marker-mid="url(#disc)" marker-end="url(#disc)" points="
      {% for row in example_results %}{{
        row.n | times: 50
      }} {{
        row.imp | times: -200 | plus: 200
      }}
      {% endfor %}" />
    <polyline stroke="none" fill="none" marker-start="url(#+)" marker-mid="url(#+)" marker-end="url(#+)" points="
      {% for row in example_results %}{{
        row.n | plus: 0.5 | times: 50
      }} {{
        row.impalt | times: -200 | plus: 200
      }}
      {% endfor %}" />
    <text x="60" y="80">
      <tspan class="katex"><tspan class="mathnormal">x</tspan></tspan><tspan dy="4" style="font-size:60%;">IMP</tspan>
    </text>
    <text x="40" y="130">
      <tspan class="katex"><tspan class="mathnormal">x</tspan></tspan><tspan dy="4" style="font-size:60%;">TR</tspan>
    </text>
    <rect x="0" width="275" y="0" height="200" stroke="#000" fill="none" />
    <g text-anchor="middle" dominant-baseline="hanging">
      <line x1="0" y1="200" x2="0" y2="205" stroke="#000" />
      <text x="0" y="207">0</text>
      <line x1="50" y1="200" x2="50" y2="205" stroke="#000" />
      <text x="50" y="207">1</text>
      <line x1="100" y1="200" x2="100" y2="205" stroke="#000" />
      <text x="100" y="207">2</text>
      <line x1="150" y1="200" x2="150" y2="205" stroke="#000" />
      <text x="150" y="207">3</text>
      <line x1="200" y1="200" x2="200" y2="205" stroke="#000" />
      <text x="200" y="207">4</text>
      <line x1="250" y1="200" x2="250" y2="205" stroke="#000" />
      <text x="250" y="207">5</text>
      <text x="138.5" y="225">
        <tspan class="katex"><tspan class="mathnormal">n</tspan></tspan>
      </text>
    </g>
    <g text-anchor="end" dominant-baseline="middle">
      <line x1="0" y1="0" x2="-5" y2="0" stroke="#000" />
      <text x="-7" y="0">1</text>
      <line x1="0" y1="100" x2="-5" y2="100" stroke="#000" />
      <text x="-7" y="100">0.5</text>
      <line x1="0" y1="200" x2="-5" y2="200" stroke="#000" />
      <text x="-7" y="200">0</text>
    </g>
  </g>
</svg>

## Application to implicit ODEs

Now consider a differential equation
```equation label=implicit_ode
g(x, \dot{x}) = 0
```
and assume there exists a function `$f$` such that 
```equation label=implicitly_defined_f
g(x, f(x)) = 0,
```
i.e. we have the corresponsing ODE `$\dot{x}=f(x)$`, but `$f$` may not be
available in closed form, so evaluating it may require expensive numerical
solving. Of course, with this definition of `$f$`,
{% eqref "implicitly_defined_f" %} should also hold for `$x^n$`, so that we have
```equation
g(x^n, f(x^n)) = 0.
```
Now we may substitute {% eqref "trapzstate" %} and {% eqref "trapzf" %} and
obtain
```equation
g\left(\tfrac12 \big(x^{n+\frac12} + x^{n-\frac12}\big),
\tfrac{1}{T}\big(x^{n+\frac12} - x^{n-\frac12}\big)\right) = 0
```
which can be used to produce a sequence `$x^{n+\frac12}$` by repeated numerical
solution and then recover `$x^n$` with {% eqref "trapzstate" %}. We have thus
applied the trapezoidel rule to {% eqref "implicit_ode" %} without requiring
`$f$`. In fact, this way we can apply the trapezoidal rule to
differential-algebraic equations (DAEs), which are of the form {% eqref
"implicit_ode" %} but where no solution function `$f$` exists. However, the
trapezoidal rule may then loose its convenient properties e.g. with regard to
stability.

For the implicit midpoint rule, we can argue similarly starting from
```equation
g(x^{n+\frac12}, f(x^{n+\frac12})) = 0
```
and inserting {% eqref "impstate" %} and
`$f\big(x^{n+\tfrac12}\big) = \tfrac{1}{T}(x^{n+1} - x^n)$` from {% eqref "impforw" %}
to obtain
```equation
g\left(\tfrac{1}{2}(x^{n+1} + x^n), \tfrac{1}{T}(x^{n+1} - x^n)\right) = 0.
```
Like for the trapezoidal rule, `$f$` is not required explicitly, but if it does
not exist, the usual properties of the implicit midpoint rule may no longer
hold.