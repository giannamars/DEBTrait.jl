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
Soil                = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))
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

CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/biocrunch/Acetate_cell_size_variation.csv", Acetate_dSize)

#
n_microbes          = 1
n_substrates        = 1
L_DNA               = 3726411*ones(1)
rRNA                = DEBTrait.genome_size_to_rRNA_copy_number(L_DNA)
min_gentime         = rand(min_gentime_distr, n_microbes)
Isolate             = BaseGenome(L_DNA, rRNA, min_gentime, 1.0*ones(n_substrates,n_microbes))
D_S                 = DEBTrait.aqueous_diffusivity(59.04402*ones(1))
Monomer             = BaseMonomer(D_S, 6*ones(n_substrates), 3*ones(n_substrates), 0.3*ones(n_substrates))
Soil                = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))
MetabolismC         = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
MetabolismC.y_EM    = ones(n_microbes)
MetabolismC.y_EX    = ones(n_microbes)
MetabolismC.α       = zeros(n_microbes)
AssimilationC       = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
#
Soil.pct_sand       = 70.0
Soil.pct_clay       = 30.0
Soil.N_cell         = 10.0*ones(n_microbes)


K_D_all = zeros(100)
s_sats = collect(range(0.01, 1.0, length=100))
for i in 1:100
    Soil.s_sat          = s_sats[i]
    AssimilationC       = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
    K_D_all[i]          = AssimilationC.K_D[1]
end

Acetate_Geobacter_dSat        = DataFrame()
Acetate_Geobacter_dSat.V_cell = fill(L_DNA[1], 100)
Acetate_Geobacter_dSat.min_gentime = fill(min_gentime[1], 100)
Acetate_Geobacter_dSat.transporter_density = fill(AssimilationC.ρ_p[1], 100)
N_SB = AssimilationC.N_SB[1,:]*100
Acetate_Geobacter_dSat.V_max = fill(N_SB[1], 100)
Acetate_Geobacter_dSat.s_sat = s_sats
Acetate_Geobacter_dSat.K = K_D_all

CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/biocrunch/Acetate_Geobacter_s_sat_variation.csv", Acetate_Geobacter_dSat)
