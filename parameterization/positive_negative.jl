using DEBTrait
using CSV, DataFrames, Statistics

n_polymers              = 0
n_monomers              = 6
n_microbes              = 2
n_enzymes               = 1
n_minerals              = 0


df_isolates             = CSV.read("/Users/glmarschmann/.julia/dev/DEBTrait/data/microbes2traits.csv", DataFrame, missingstring = "N/A")
df_monomers             = CSV.read("/Users/glmarschmann/.julia/dev/DEBTrait/files/thermodynamic_props_class.csv", DataFrame, missingstring = "N/A")
df_p                    = filter(row -> row[:Rhizosphere_response] == "positive", df_isolates)
df_n                    = filter(row -> row[:Rhizosphere_response] == "negative", df_isolates)

L_DNA                   = [median(df_p.Genome_size), median(df_n.Genome_size)]
rRNA                    = [median(df_p.rRNA_genes), median(df_n.rRNA_genes)]
min_gentime             = [median(df_p.Min_gen_time), median(df_n.Min_gen_time)]

z_sugars_p              = median(convert(Array{Float64,1}, df_p.z_sugars))*ones(1,1)
z_organics_p            = median(convert(Array{Float64,1}, df_p.z_organic_acids))*ones(1,1)
z_aminos_p              = median(convert(Array{Float64,1}, df_p.z_amino_acids))*ones(1,1)
z_fattys_p              = median(convert(Array{Float64,1}, df_p.z_fatty_acids))*ones(1,1)
z_nucleos_p             = median(convert(Array{Float64,1}, df_p.z_nucleotides))*ones(1,1)
z_auxins_p              = median(convert(Array{Float64,1}, df_p.z_auxins))*ones(1,1)
genome_distr_p          = vcat(z_sugars_p, z_organics_p, z_aminos_p, z_fattys_p, z_nucleos_p, z_auxins_p)

z_sugars_n              = median(convert(Array{Float64,1}, df_n.z_sugars))*ones(1,1)
z_organics_n            = median(convert(Array{Float64,1}, df_n.z_organic_acids))*ones(1,1)
z_aminos_n              = median(convert(Array{Float64,1}, df_n.z_amino_acids))*ones(1,1)
z_fattys_n              = median(convert(Array{Float64,1}, df_n.z_fatty_acids))*ones(1,1)
z_nucleos_n             = median(convert(Array{Float64,1}, df_n.z_nucleotides))*ones(1,1)
z_auxins_n              = median(convert(Array{Float64,1}, df_n.z_auxins))*ones(1,1)
genome_distr_n          = vcat(z_sugars_n, z_organics_n, z_aminos_n, z_fattys_n, z_nucleos_n, z_auxins_n)

genome_distr            = hcat(genome_distr_p, genome_distr_n)


Isolate                 = BaseGenome(L_DNA, rRNA, min_gentime, genome_distr)
Monomer                 = BaseMonomer(df_monomers.D_S, df_monomers.N_C, df_monomers.N_N, df_monomers.y_DE)
Soil                    = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))


MetabolismC             = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
MetabolismC.y_EM        = ones(n_microbes)
MetabolismC.y_EX        = ones(n_microbes)
MetabolismC.α           = zeros(n_microbes)

AssimilationC           = PAssimilation(Isolate, Monomer, Soil, MetabolismC)

Setup                   = SetupBatchC(Isolate, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
p                       = ParamsBatchC(Setup, Monomer, MetabolismC, AssimilationC)
