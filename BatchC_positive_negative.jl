using DEBTrait
using CSV, DataFrames, Statistics
using DifferentialEquations, Plots

include("/Users/glmarschmann/.julia/dev/DEBTrait/parameterization/positive_negative.jl")

u0 = zeros(p.dim)
u0[1+p.n_polymers:p.n_polymers+p.n_monomers] = 25.0*ones(p.n_monomers)
u0[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes] = p.Setup.u0_reserves
u0[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes] = p.Setup.u0_structures

p.PMetabolismC.α = [0.01, 0.3]
p.PMetabolismC.y_EX = [1.4, 2.8]

tspan = (0.0,3.3)
prob  = ODEProblem(rhs_BatchC!,u0,tspan,p)
sol   = solve(prob, alg_hints=[:stiff])


plot(sol, vars=1+p.n_polymers:p.n_polymers+p.n_monomers)
plot(sol, vars=1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes)
plot(sol, vars=1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes)
plot(sol, vars=1+p.n_polymers+p.n_monomers+2*p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes)
plot(sol, vars=1+p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes+p.n_microbes)



BGE_pos = [DEBTrait.calc_BGE(zeros(p.dim), sol.u[i], p, 0.0)[1] for i in 1:size(sol.t,1)]
BGE_neg = [DEBTrait.calc_BGE(zeros(p.dim), sol.u[i], p, 0.0)[2] for i in 1:size(sol.t,1)]
BGE_median = [median(BGE_pos), median(BGE_neg)]
plot(sol.t, BGE_pos)
plot!(sol.t, BGE_neg)



r_pos = [DEBTrait.calc_growth_rate(sol.u[i], p)[1] for i in 1:size(sol.t,1)]
r_neg = [DEBTrait.calc_growth_rate(sol.u[i], p)[2] for i in 1:size(sol.t,1)]
r_median = [median(r_pos), median(r_neg)]
plot(sol.t, r_pos)
plot!(sol.t, r_neg)

scatter(r_median, BGE_median)


x_pos = [DEBTrait.calc_growth_production(sol.u[i], p)[1][1] for i in 1:size(sol.t,1)]
x_neg = [DEBTrait.calc_growth_production(sol.u[i], p)[1][2] for i in 1:size(sol.t,1)]
plot(sol.t, x_pos)
plot!(sol.t, x_neg)
