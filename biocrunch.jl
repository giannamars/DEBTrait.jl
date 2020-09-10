using DEBTrait
using Distributions
using CSV, DataFrames
using Random
Random.seed!(123)

n_microbes          = 100
n_substrates        = 1
#
L_DNA               = collect(range(1e6, 8e6, length=n_microbes))
rRNA                = DEBTrait.genome_size_to_rRNA_copy_number(L_DNA)
min_gentime_distr   = Gamma(2.486, 0.97)
min_gentime         = rand(min_gentime_distr, n_microbes)
Isolate             = BaseGenome(L_DNA, rRNA, min_gentime, 1.0*ones(n_substrates,n_microbes))
D_S                 = DEBTrait.aqueous_diffusivity(59.04402*ones(1))
Monomer             = BaseMonomer(D_S[1]*ones(n_substrates), 6*ones(n_substrates), 3*ones(n_substrates), 0.3*ones(n_substrates))
Soil                = BaseSoil(0.0, 0.0, 1.0, 10*ones(n_microbes))
MetabolismC         = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
MetabolismC.y_EM    = ones(n_microbes)
MetabolismC.y_EX    = ones(n_microbes)
MetabolismC.α       = zeros(n_microbes)
AssimilationC       = PAssimilation(Isolate, Monomer, Soil, MetabolismC)

Acetate_dSize        = DataFrame()
Acetate_dSize.V_cell = DEBTrait.genome_size_to_cell_volume(L_DNA)
Acetate_dSize.min_gentime = min_gentime
Acetate_dSize.transporter_density = AssimilationC.ρ_p[1,:]
Acetate_dSize.K      = AssimilationC.K_D_0[1,:]
Acetate_dSize.V_max = AssimilationC.N_SB[1,:]*100

CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/Acetate_cell_size_variation.csv", Acetate_dSize)

#
n_microbes          = 1
n_substrates        = 1
L_DNA               = 3726411*ones(1)
rRNA                = DEBTrait.genome_size_to_rRNA_copy_number(L_DNA)
min_gentime         = rand(min_gentime_distr, n_microbes)
Isolate             = BaseGenome(L_DNA, rRNA, min_gentime, 1.0*ones(n_substrates,n_microbes))
D_S                 = DEBTrait.aqueous_diffusivity(59.04402*ones(1))
Monomer             = BaseMonomer(D_S, 6*ones(n_substrates), 3*ones(n_substrates), 0.3*ones(n_substrates))
Soil                = BaseSoil(0.0, 0.0, 1.0, 10*ones(n_microbes))
MetabolismC         = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
MetabolismC.y_EM    = ones(n_microbes)
MetabolismC.y_EX    = ones(n_microbes)
MetabolismC.α       = zeros(n_microbes)
AssimilationC       = PAssimilation(Isolate, Monomer, Soil, MetabolismC)

Soil.pct_sand    = 30.0
Soil.pct_clay    = 70.0
Soil.s_sat       = 0.01
Soil.N_cell      = 100.0*ones(n_microbes)
AssimilationC       = PAssimilation(Isolate, Monomer, Soil, MetabolismC)


AssimilationC.K_D
