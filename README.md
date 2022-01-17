# Fast Sweeping & Fast Marching methods for the solution of eikonal equations

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://triscale-innov.github.io/FastSweeping.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://triscale-innov.github.io/FastSweeping.jl/dev)
[![Build Status](https://github.com/triscale-innov/FastSweeping.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/triscale-innov/FastSweeping.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/triscale-innov/FastSweeping.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/triscale-innov/FastSweeping.jl)

Julia implementations of solvers for general Eikonal equations of the form

$$\begin{align}
\left\vert\nabla t\right\vert^2 = \frac{1}{c(x,y)}, & \forall (x,y)\in\Omega\subset\mathbb{R}^2,\\
t(x_0,y_0) = 0, & \forall (x_0,y_0)\in\Gamma\subset\Omega,
\end{align}$$

where $\Omega$ is a rectangular, 2D spatial domain, and $t(x,y)$ represents the
first arrival time at point $(x,y)$ of a front moving at speed $c(x,y)$ and
originating from $\Gamma$.

This package provides implementations for two methods

- Fast Sweeping Method (FSM)  [1]
- Fast Marching Method (FMM)  [2]


[1] Zhao, Hongkai (2005-01-01). "A fast sweeping method for Eikonal equations". Mathematics of Computation. 74 (250): 603–627. [DOI: 10.1090/S0025-5718-04-01678-3](https://doi.org/10.1090%2FS0025-5718-04-01678-3)<br/>
[2] J.A. Sethian. A Fast Marching Level Set Method for Monotonically Advancing Fronts, Proc. Natl. Acad. Sci., 93, 4, pp.1591--1595, 1996. [PDF](https://math.berkeley.edu/~sethian/2006/Papers/sethian.fastmarching.pdf)
