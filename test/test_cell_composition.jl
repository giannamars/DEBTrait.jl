using DEBTrait, Test
using CSV, DataFrames, Statistics
using Plots

df_isolates      = CSV.read("/Users/glmarschmann/.julia/dev/DEBTrait/data/microbes2traits.csv", DataFrame, missingstring = "N/A")
df_isolates      = filter(row -> row[:Rhizosphere_response] == "negative" || row[:Rhizosphere_response] == "positive" , df_isolates)

L_DNA            = convert(Array{Float64,1}, df_isolates.Genome_size)
rRNA             = DEBTrait.genome_size_to_rRNA_copy_number(L_DNA)       # need to cap this at 20
rRNA_iso         = df_isolates.rRNA_genes
min_gentime      = df_isolates.Min_gen_time
Isolate          = BaseGenome(L_DNA, rRNA, min_gentime, 1.0*ones(1,27), 1.0*ones(27))
Cell_Composition = PCellComposition{BaseGenome}(Isolate)

rRNA_protein_ratio = Cell_Composition.V_r./Cell_Composition.V_p
y_EV             = DEBTrait.reserve_yield_constant(Cell_Composition)
