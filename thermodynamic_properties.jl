using DEBTrait
using CSV, DataFrames, Statistics

## this is all still very clumsy

df_monomers = CSV.read("/Users/glmarschmann/.julia/dev/DEBTrait/data/monomers_all.csv", DataFrame, missingstring="N/A")
chemFormBiom = [1, 1.8, 0.2, 0.5, 0, 0, 0]

n_monomers = size(df_monomers,1)
dGcox = zeros(n_monomers)
for i in 1:n_monomers
    elementstring = df_monomers.Formula[i]
    out = DEBTrait.ThermoStoichWizard.get_lambda(elementstring, chemFormBiom)
    dGcox[i] = out[2][3]
end
df_monomers.dGcox = dGcox

dGcat = zeros(n_monomers)
for i in 1:n_monomers
    elementstring = df_monomers.Formula[i]
    out = DEBTrait.ThermoStoichWizard.get_lambda(elementstring, chemFormBiom)
    dGcat[i] = out[2][5]
end
df_monomers.dGcat = dGcat

dGAn = zeros(n_monomers)
for i in 1:n_monomers
    elementstring = df_monomers.Formula[i]
    out = DEBTrait.ThermoStoichWizard.get_lambda(elementstring, chemFormBiom)
    dGAn[i] = out[2][8]
end
df_monomers.dGAn = dGAn

λ = zeros(n_monomers)
for i in 1:n_monomers
    elementstring = df_monomers.Formula[i]
    out = DEBTrait.ThermoStoichWizard.get_lambda(elementstring, chemFormBiom)
    λ[i] = out[1][1]
end
df_monomers.lambda = λ

η    = 0.43
y_DE = @. (λ*η*dGcat)/dGcox
df_monomers.y_DE = y_DE

N_C = zeros(n_monomers)
for i in 1:n_monomers
    elementstring = df_monomers.Formula[i]
    N_C[i]        = DEBTrait.ThermoStoichWizard.extract_composition(elementstring)[1]
end
df_monomers.N_C = N_C

N_N = zeros(n_monomers)
for i in 1:n_monomers
    elementstring = df_monomers.Formula[i]
    N_N[i]        = DEBTrait.ThermoStoichWizard.extract_composition(elementstring)[3]
end
df_monomers.N_N = N_N

D_S           = DEBTrait.aqueous_diffusivity(convert(Array{Float64,1}, df_monomers.Molecular_weight))

df            = DataFrame()
df.y_DE       = y_DE
df.D_S        = D_S
df.N_C        = N_C
df.N_N        = N_N
df.Ontology   = df_monomers.Ontology
CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/thermodynamic_props_all.csv", df)

############ calculate median values for monomer classes
tmp           = convert(Array{Int64,1},1:n_monomers)
sugars        = tmp[df_monomers.Ontology .== "sugar"]
organics      = tmp[df_monomers.Ontology .== "organic acid"]
aminos        = tmp[df_monomers.Ontology .== "amino acid"]
fattys        = tmp[df_monomers.Ontology .== "fatty acid"]
nucleos       = tmp[df_monomers.Ontology .== "nucleotide"]
auxins        = tmp[df_monomers.Ontology .== "auxin"]

y_DE          = vcat(median(y_DE[sugars]),median(y_DE[organics]), median(y_DE[aminos]), median(y_DE[fattys]), median(y_DE[nucleos]), median(y_DE[auxins]))
D_S           = vcat(median(D_S[sugars]),median(D_S[organics]), median(D_S[aminos]), median(D_S[fattys]), median(D_S[nucleos]), median(D_S[auxins]))
N_C           = vcat(median(N_C[sugars]),median(N_C[organics]), median(N_C[aminos]), median(N_C[fattys]), median(N_C[nucleos]), median(N_C[auxins]))
N_N           = vcat(median(N_N[sugars]),median(N_N[organics]), median(N_N[aminos]), median(N_N[fattys]), median(N_N[nucleos]), median(N_N[auxins]))

df_classes    = DataFrame()
df_classes.y_DE = y_DE
df_classes.D_S = D_S
df_classes.N_C = N_C
df_classes.N_N = N_N
df_classes.Ontology = ["sugar", "organic acid", "amino acid", "fatty acid", "nucletotide", "auxin"]

CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/thermodynamic_props_class.csv", df_classes)
