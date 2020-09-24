n_microbes      = 1
n_reserves      = 2
E               = fill(0.5,n_microbes*n_reserves)
E               = reshape(E, n_reserves, n_microbes)
V               = fill(1.0, n_microbes)

Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes), 1.0*ones(n_microbes), 1.0*ones(1,n_microbes), 1.0*ones(n_microbes))

p               = PMetabolismCN{BaseGenome, Array{Float64,2}, Array{Float64,1}}(Isolate)
p.α             = fill(0.1, n_microbes)
p.α_n           = fill(0.1, n_microbes)
p.y_EM          = ones(n_reserves, n_microbes)
p.y_EX          = ones(n_reserves, n_microbes)
p.y_EX[1,1]     = 0.1
p.y_EX[2,1]     = 0.3
r               = growth!(zeros(size(E,2)), E, V, p)

x, rM_CO2, rEG_rej, rEX_rej = growth_production(r, E, V, p)
