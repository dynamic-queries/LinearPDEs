include("includes.jl")
pyplot()

function Laplacian(x)
    dx = x[2]-x[1]
    k = -2*ones(ComplexF64,length(x))
    k1 = ones(ComplexF64,length(x)-1)
    D = diagm(0=>k) + diagm(1=>k1) + diagm(-1=>k1)
    D = (1/dx)^2 .* D
    D
end 

function initialize(ψ1, ψ2)
    (ψ1 + ψ2)/sqrt(2)
end 

# Problem parameters
a = 1e-9 # m
N = 100
z = LinRange(-a,a,N)
z0 = a/(4*sqrt(2))
Vb = 1*1.6e-19 # J
# E = 1e9 # V/m
E = 0
ω = 2π*100*1e12 # 1/s
h = 1.38e-34 # JK
me = 9.1e-31 # kg

# Define potential
V(x) = Vb * (-0.25*(x/z0)^2 + 0.015625*(x/z0)^4)

# Eval potential
U = V.(z)


# Assemble hamiltonian
H = (-h^2/(2*me))*Laplacian(z) + diagm(U)
E, ϕ = eigen(H)
ψ₀ = initialize(ϕ[:,1],ϕ[:,2])

plot(z/1e-9,U/(8*Vb),title="Double well potential",label=false,xlabel=L"z [nm]",ylabel=L"U/(8V_b) \:[no\:unit]")
plot!(twinx(), z/1e-9,abs.(ψ₀),label=false, c=:red, ylabel=L"|\psi| \:[no\:unit]",ylim=[0,1])
savefig("figures/schrodinger_initial.png")

# Forcing function
function SE!(du, u, p, t)
    Vt = 1.6e-19 * E*sin(ω*t) .* ones(length(u)) 
    du .= (-1im/h).*(H - diagm(0=>Vt))*u
    nothing
end 

# Setup time integration
params = []
tspan = (0, 1e-13)
δt = 1e-15
t = LinRange(tspan[1],tspan[2],N)
init = ψ₀
prob = ODEProblem(SE!, init, tspan, params)

sol = solve(prob, Tsit5(),dt=δt, saveat=t)
solution = Array(sol)
anim = @animate for i=1:size(solution,2)
    plot(z/1e-9,U/(8*Vb),title="Time step = $(i)",label=false,xlabel=L"z [nm]",ylabel=L"U/(8V_b) \:[no\:unit]")
    ψ₀ = solution[:,i]
    plot!(twinx(), z/1e-9,abs.(ψ₀),label=false, c=:red, ylabel=L"|\psi| \:[no\:unit]",ylim=[0,1])
end 
gif(anim, "figures/SE_explicit.gif", fps=5)

sol = solve(prob, Trapezoid(autodiff=false),dt=δt, saveat=t)
solution = Array(sol)
anim = @animate for i=1:size(solution,2)
    plot(z/1e-9,U/(8*Vb),title="Time step = $(i)",label=false,xlabel=L"z [nm]",ylabel=L"U/(8V_b) \:[no\:unit]")
    ψ₀ = solution[:,i]
    plot!(twinx(), z/1e-9,abs.(ψ₀),label=false, c=:red, ylabel=L"|\psi| \:[no\:unit]",ylim=[0,1])
end 
gif(anim, "figures/SE_implicit.gif", fps=5)

# Figure
gr()
f1 = Plots.heatmap(z,t,abs.(solution),c=:summer,xlabel=L"x",ylabel=L"t")
title!("Schrodinger equation with Double well potential")
savefig("figures/contour_SE.png")
savefig("figures/contour_SE.svg")

# Write data
file = h5open("data/SE","w")
file["README"] = "Solution to 1D Schroginger equation with double well potential. The system is known to be stiff."
file["x"] = Vector(z)
file["t"] = Vector(t)
file["ψxt"] = solution
close(file)