using DEBTrait
using CSV, DataFrames, Statistics

df_isolates             = CSV.read("/Users/glmarschmann/.julia/dev/DEBTrait/data/microbes2traits.csv", DataFrame, missingstring = "N/A")
df_isolates             = filter(row -> row[:Rhizosphere_response] == "positive" || row[:Rhizosphere_response] == "negative", df_isolates)


L_DNA                   = df_isolates.Genome_size
rRNA                    = df_isolates.rRNA_genes
min_gentime             = df_isolates.Min_gen_time
hydrolases              = df_isolates.hydrolases
#
z_sugars_p              = reshape(df_isolates.z_sugars, 1, 27)
z_organics_p            = reshape(df_isolates.z_organic_acids, 1, 27)
z_aminos_p              = reshape(df_isolates.z_amino_acids, 1, 27)
z_fattys_p              = reshape(df_isolates.z_fatty_acids, 1, 27)
z_nucleos_p             = reshape(df_isolates.z_nucleotides, 1, 27)
z_auxins_p              = reshape(df_isolates.z_auxins, 1, 27)
genome_distr_p          = vcat(z_sugars_p, z_organics_p, z_aminos_p, z_fattys_p, z_nucleos_p, z_auxins_p)

Isolate                 = BaseGenome(L_DNA, rRNA, min_gentime, genome_distr_p[1]*ones(n_monomers,n_microbes), hydrolases)

MetabolismC             = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)

isolate_id              = collect(1:27)
id_positive             = isolate_id[df_isolates.Rhizosphere_response .== "positive"]
id_negative             = isolate_id[df_isolates.Rhizosphere_response .== "negative"]

median(MetabolismC.y_EV[id_positive])
median(MetabolismC.y_EV[id_negative])
