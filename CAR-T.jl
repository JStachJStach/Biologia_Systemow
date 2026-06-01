## środowisko
using Pkg
Pkg.activate("env/")
Pkg.instantiate()
using CairoMakie
import DifferentialEquations as DE
using DataFrames
using LinearAlgebra
using Random
##parametry
u0 = [10; 10]
tspan = (0.0, 13.0)
parameters = (100, 0.1, 0.05, 0.01, 0.5)


## zestaw równań 
function f!(du, u, p, t)
    T, E = u
    kin, gammaT, kkill, P, gammaE = p#100, 0.1, 0.05, 0.01, 0.5

    #if t > 10
    #    kkill = 0.5
    #end
    du[1] = dTdt = kin - gammaT * T - kkill * E * T
    du[2] = dEdt = P * E * T - gammaE * E
end




## rozwiązanie i plotting
prob = DE.ODEProblem(f!, u0, tspan, parameters)
sol = DE.solve(prob, saveat=0.01)
Ts = [u[1] for u in sol.u]
Es = [u[2] for u in sol.u]



fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, sol.t, Ts, label="T")
lines!(ax, sol.t, Es, label="E")

display(fig)


#############################
# rozwiązanie stochastyczne #
#############################
##warunki początkowe:
kin, gammaT, kkill, P, gammaE = parameters
T, E = u0

T_vec = []
E_vec = []
t_vec = []
last_t = 0
## symulacja:
for i in 1:3000
    r1 = rand(Float64)
    r2 = rand(Float64)
    v1 = kin                 #synteza T
    v2 = P * E * T           #synteza E
    v3 = gammaT * T + kkill * E * T #degradajca T
    v4 = gammaE * E            #degradajca E


    h = v1 + v2 + v3 + v4
    p1 = v1 / h
    p2 = v2 / h
    p3 = v3 / h
    p4 = v4 / h
    last_t += -log(r1) / h

    #synteza
    if r2 < p1
        T += 1
    end
    if r2 < p2
        E += 1
    end
    #degradacja
    if r2 < p3
        T -= 1
    end
    if r2 < p4
        E -= 1
    end
    append!(T_vec, T)
    append!(E_vec, E)
    append!(t_vec, last_t)
end


lines!(ax, t_vec, T_vec, label="T_stoch")
lines!(ax, t_vec, E_vec, label="E_stoch")
axislegend(ax, position=:rb)
display(fig)

