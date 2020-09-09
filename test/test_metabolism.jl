using DEBTrait, Test
using Random
Random.seed!(123)

################# (1 microbe, 1 reserve), no enzyme production
n_microbes      = 1
E               = rand(n_microbes)
V               = rand(n_microbes)
Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes))
p               = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
p.y_EM          = fill(1.0, n_microbes)
p.y_EX          = fill(1.0, n_microbes)
r               = growth!(zeros(size(E,1)), E, V, p)
r_analytic      = @. (p.k_E*E/V - p.k_M*p.y_EM)/(E/V + p.y_EV)
@test r ≈ r_analytic atol=1e-10

out             = growth_production(r, E, V, p)
x_analytic      = @. p.α*p.y_EV*r_analytic/p.y_EX
@test out[1] ≈ x_analytic atol=1e-10

################# (1 microbe, 2 reserves)
n_reserves      = 2
E               = fill(0.5,n_microbes*n_reserves)
E               = reshape(E, n_reserves, n_microbes)
V               = fill(1.0, n_microbes)
p               = PMetabolismCN{BaseGenome, Array{Float64,2}, Array{Float64,1}}(Isolate)
p.y_EX          = fill(1.0, n_reserves, n_microbes)
p.k_M           = fill(1.0, n_microbes)
p.α             = fill(0.0, n_microbes)
p.y_EM[:,1]     = [0.3, 1.0]
p.y_EV[:,1]     = [1.2, 1.0]
p.k_E           = [1.0]
r               = growth!(zeros(size(E,2)), E, V, p)
@test r[1] ≈ -0.333 atol=1e-2
out             = growth_production(r, E, V, p)
@test out[5][:,1] ≈ [0.3667, 0.0] atol=1e-4

################# (2 microbes, 2 reserves)
n_microbes      = 2
n_reserves      = 2
Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes))
E               = fill(0.5,n_microbes*n_reserves)
E               = reshape(E, n_reserves, n_microbes)
V               = fill(1.0, n_microbes)
p               = PMetabolismCN{BaseGenome, Array{Float64,2}, Array{Float64,1}}(Isolate)
r               = growth!(zeros(size(E,2)), E, V, p)
@test size(r,1) == 2
out             = growth_production(r, E, V, p)
@test size(out[1]) == (n_reserves,n_microbes)
