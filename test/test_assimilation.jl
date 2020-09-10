using DEBTrait, Test
using Random
Random.seed!(123)

######## benchmark K_D against Tang et al. (2019)
D_S     = 1.4e-9*ones(1)
N_p     = 3000*ones(1,1)
r_c     = 1e-6
V_cell  = (4/3)*pi*r_c^3*ones(1)
ρ_p     = DEBTrait.transporter_density(V_cell, N_p)
K_D_0   = DEBTrait.specific_reference_affinity(ρ_p, V_cell, D_S)
@test K_D_0[1] ≈ 5.7969874109e-05 atol=1e-6
#
pct_sand = 70.0
pct_clay = 30.0
s_sat    = 0.10101010101010102
N_cell   = 10.0*ones(1)
Soil = BaseSoil(pct_sand, pct_clay, s_sat, N_cell)
K_D     = DEBTrait.affinity_constant(ρ_p, V_cell, D_S, Soil.pct_sand, Soil.pct_clay, Soil.s_sat, Soil.N_cell)
@test K_D[1] ≈ 1.447412099595345 atol=1e-6
#
Soil.s_sat = 0.8080808080808082
K_D     = DEBTrait.affinity_constant(ρ_p, V_cell, D_S, Soil.pct_sand, Soil.pct_clay, Soil.s_sat, Soil.N_cell)
@test K_D[1] ≈ 9.390492906747054e-05 atol=1e-6
#
N_SB    = DEBTrait.monomer_uptake_sites(V_cell, ρ_p)

########
n_microbes      = 1
n_substrates    = 1
Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes), 1.0*ones(n_microbes), 1.0*ones(n_substrates,n_microbes))
Monomer         = BaseMonomer(1e-10*ones(n_substrates), 6*ones(n_substrates), 3*ones(n_substrates), 0.3*ones(n_substrates))
Soil            = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))
MetabolismC     = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
MetabolismC.y_EM = ones(n_microbes)
MetabolismC.y_EX = ones(n_microbes)
MetabolismC.α    = zeros(n_microbes)
AssimilationC   = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
@test size(AssimilationC.ρ_p) == (n_substrates, n_microbes)
@test size(AssimilationC.N_SB) == (n_substrates, n_microbes)
@test size(AssimilationC.K_D) == (n_substrates, n_microbes)
@test AssimilationC.K_D_0 ≈ AssimilationC.K_D atol = 1e-4
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
