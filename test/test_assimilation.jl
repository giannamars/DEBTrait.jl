using DEBTrait, Test
using Random
Random.seed!(123)

######## benchmark K_D against Tang et al. (2019)
D_S     = 5e-10*ones(1)
N_p     = 3140*ones(1,1)
r_c     = 1e-6
V_cell  = (4/3)*pi*r_c^3*ones(1)
ρ_p     = DEBTrait.transporter_density(V_cell, N_p)
K_D     = DEBTrait.specific_reference_affinity(ρ_p, V_cell, D_S)
@test K_D ≈ 1.66e-4*ones(1,1) atol=1e-7
N_SB    = DEBTrait.monomer_uptake_sites(V_cell, ρ_p)

########
n_microbes      = 1
n_substrates    = 1
Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes), 1.0*ones(n_microbes), 1.0*ones(n_substrates,n_microbes))
Monomer         = BaseMonomer(1e-10*ones(n_substrates), 6*ones(n_substrates), 3*ones(n_substrates), 0.3*ones(n_substrates))
Soil            = BaseSoil(0.0, 0.0, 1.0, 10*ones(n_microbes))
MetabolismC     = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
MetabolismC.y_EM = ones(n_microbes)
MetabolismC.y_EX = ones(n_microbes)
MetabolismC.α    = zeros(n_microbes)
AssimilationC   = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
@test size(AssimilationC.ρ_p) == (n_substrates, n_microbes)
@test size(AssimilationC.N_SB) == (n_substrates, n_microbes)
@test size(AssimilationC.K_D) == (n_substrates, n_microbes)
@test AssimilationC.K_D_0 ≈ AssimilationC.K_D atol = 1e-8
Soil.s_sat = 0.05
AssimilationC   = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
@test AssimilationC.K_D_0[1] < AssimilationC.K_D[1]


n_microbes      = rand(1:10)
n_substrates    = rand(1:10)
Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes), 1.0*ones(n_microbes), 1.0*ones(n_substrates,n_microbes))
Monomer         = BaseMonomer(1e-10*ones(n_substrates), 6*ones(n_substrates), 3*ones(n_substrates), 0.3*ones(n_substrates))
Soil            = BaseSoil(0.0, 0.0, 1.0, 10*ones(n_microbes))
MetabolismC     = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
MetabolismC.y_EM = ones(n_microbes)
MetabolismC.y_EX = ones(n_microbes)
MetabolismC.α    = zeros(n_microbes)
AssimilationC   = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
@test size(AssimilationC.ρ_p) == (n_substrates, n_microbes)
@test size(AssimilationC.N_SB) == (n_substrates, n_microbes)
@test size(AssimilationC.K_D) == (n_substrates, n_microbes)
