include("includes.jl")

Nx = 400
Nt = 400
x = Vector(LinRange(0.0, π, Nx))
t = Vector(LinRange(0.0, 10.0, Nt))
dx = x[2]-x[1]
cx = 1.5

function initialize(x::AbstractVector)
    σ = 0.5
    μ1 = π/3
    μ2 = 2*π/3
    f(x) = exp(-((x-μ1)/σ)^2) + exp(-((x-μ2)/σ)^2)
    g(x) = -(2*((x-μ1))/σ) * exp(-((x-μ1)/σ)^2) -(2*((x-μ2))/σ) * exp(-((x-μ2)/σ)^2)
    vcat(f.(x), g.(x))
end 

function wave!(du, u, p, t)
    k = p[1]
    n = floor(Int, length(u)/2)
    for i=2:n-1
        du[n+i] = k*(u[i+1]+u[i-1]-2*u[i])
        du[i] = u[n+i]
    end 
    nothing
end

init = initialize(x)
tspan = (t[1],t[end])
params = [(cx/dx)^2]
prob = ODEProblem(wave!, init, tspan, params)
sol = solve(prob, Trapezoid(), saveat=t)
solution = Array(sol)

# Validate
anim = @animate for i=1:size(solution,2)
    plot(x, solution[1:length(x),i],title="Time step:$(i)",label=false, xlim=[0,π], ylim=[-2,2],xlabel=L"x", ylabel= L"u")
end 
gif(anim,"figures/wave.gif",fps=5)

# Figure
f1 = Plots.heatmap(x,t,solution[1:length(x),:],c=:summer,xlabel=L"x",ylabel=L"t")
title!("1D Wave equation")
savefig("figures/contour_wave.svg")

# Write data
file = h5open("data/wave","w")
file["README"] = "Solution to 1D wave equation with no source. Homogenous dirichlet boundary conditions. Initial condition with Gaussian bumps (SD=0.5) centered at means=(π/3,2π/3)"
file["x"] = x
file["t"] = t
file["uxt"] = solution
close(file)