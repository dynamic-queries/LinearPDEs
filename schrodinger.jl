using OrdinaryDiffEq
using Plots
using Measures

a = 1e-9 # m
N = 50
z = LinRange(-a,a,N)
z0 = a/(4*sqrt(2))
Vb = 1*1.6e-19 # J

# Define potential
V(x) = Vb * (-0.25*(x/z0)^2 + 0.015625*(x/z0)^4)

# Eval potential
U = V.(z)

fig = plot()