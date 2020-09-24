using DEBTrait
using CSV, DataFrames, Statistics
using DifferentialEquations, Plots

n_polymers              = 0
n_monomers              = 6
n_microbes              = 1
n_enzymes               = 1
n_minerals              = 0

df_isolates             = CSV.read("/Users/glmarschmann/.julia/dev/DEBTrait/data/microbes2traits.csv", DataFrame, missingstring = "N/A")
df_isolates             = filter(row -> row[:Rhizosphere_response] == "positive" || row[:Rhizosphere_response] == "negative", df_isolates)
df_monomers             = CSV.read("/Users/glmarschmann/.julia/dev/DEBTrait/files/thermodynamic_props_class.csv",DataFrame, missingstring = "N/A")


isolate_id              = collect(1:27)
monomer_id              = collect(1:6)
monomer_concentration   = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1.0]
#
BGE                     = zeros(27, 6, 13)
r                       = zeros(27, 6, 13)
m_E                     = zeros(27, 6, 13)
x                       = zeros(27, 6, 13)
#
function condition(u,t,integrator) # Event when event_f(u,t) == 0
  u[1] < 1e-8
end
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition,affect!)

i = 1
id_isolate              = isolate_id[i]
#
L_DNA                   = df_isolates.Genome_size[id_isolate]*ones(n_microbes)
rRNA                    = df_isolates.rRNA_genes[id_isolate]*ones(n_microbes)
min_gentime             = df_isolates.Min_gen_time[id_isolate]*ones(n_microbes)
#
z_sugars_p              = df_isolates.z_sugars[id_isolate]*ones(1,1)
z_organics_p            = df_isolates.z_organic_acids[id_isolate]*ones(1,1)
z_aminos_p              = df_isolates.z_amino_acids[id_isolate]*ones(1,1)
z_fattys_p              = df_isolates.z_fatty_acids[id_isolate]*ones(1,1)
z_nucleos_p             = df_isolates.z_nucleotides[id_isolate]*ones(1,1)
z_auxins_p              = df_isolates.z_auxins[id_isolate]*ones(1,1)
genome_distr_p          = vcat(z_sugars_p, z_organics_p, z_aminos_p, z_fattys_p, z_nucleos_p, z_auxins_p)
    #
j = 6
n_monomers              = 1
id_monomer              = monomer_id[j]
Isolate                 = BaseGenome(L_DNA, rRNA, min_gentime, genome_distr_p[id_monomer]*ones(n_monomers,n_microbes))
Monomer                 = BaseMonomer(df_monomers.D_S[id_monomer]*ones(n_monomers), df_monomers.N_C[id_monomer]*ones(n_monomers), df_monomers.N_N[id_monomer]*ones(n_monomers), df_monomers.y_DE[id_monomer]*ones(n_monomers))
Soil                    = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))
#
MetabolismC             = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
MetabolismC.y_EM        = ones(n_microbes)
MetabolismC.y_EX        = ones(n_microbes).*1.4
MetabolismC.α           = zeros(n_microbes)

if df_p.Rhizosphere_response[1] .== "positive"
    MetabolismC.α .= 0.01
else
    MetabolismC.α .= 0.3
end
#

AssimilationC           = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
#
Setup                   = SetupBatchC(Isolate, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
p                       = ParamsBatchC(Setup, Monomer, MetabolismC, AssimilationC)
#########################
p.n_polymers              = 0
p.n_monomers              = 1
p.n_microbes              = 1
p.n_enzymes               = 1
p.n_minerals              = 0
#
u0 = zeros(p.dim)
        for k in 1:13
            init_D0 = monomer_concentration[k]
            u0[1+p.n_polymers:p.n_polymers+p.n_monomers] = init_D0*ones(p.n_monomers)
            u0[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes] = p.Setup.u0_reserves*1e-6
            u0[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes] = p.Setup.u0_structures*1e-6
            #
            tspan = (0.0,  150.0)
            prob  = ODEProblem(rhs_BatchC!, u0, tspan, p)
            sol   = solve(prob, alg_hints=[:stiff], callback=cb)
            #
            BGE_tmp = [DEBTrait.calc_BGE(zeros(p.dim), sol.u[i], p, 0.0)[1] for i in 1:size(sol.t,1)]
            BGE[i,j,k] = median(BGE_tmp[BGE_tmp.>1e-8])
            r_tmp   = [DEBTrait.calc_growth_rate(sol.u[i], p)[1] for i in 1:size(sol.t,1)]
            r[i,j,k] = median(r_tmp)
            E = sol.u[end][1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
            V = sol.u[end][1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
            m_E_tmp = E./V
            m_E_tmp[m_E_tmp.>=1.0] .= 0.0
            m_E[i,j,k] = m_E_tmp[1]
            x_tmp = [DEBTrait.calc_growth_production(sol.u[i], p)[1][1] for i in 1:size(sol.t,1)]
            x[i,j,k] = median(x_tmp)
        end


scatter(monomer_concentration, BGE[1,1,:], xscale=:log10)
scatter!(monomer_concentration, BGE[2,1,:], xscale=:log10)

scatter(monomer_concentration, BGE[1,2,:], xscale=:log10)
scatter!(monomer_concentration, BGE[2,2,:], xscale=:log10)

scatter(monomer_concentration, BGE[1,3,:], xscale=:log10)
scatter!(monomer_concentration, BGE[2,3,:], xscale=:log10)

scatter(monomer_concentration, BGE[1,4,:], xscale=:log10)
scatter!(monomer_concentration, BGE[2,4,:], xscale=:log10)

scatter(monomer_concentration, BGE[1,5,:], xscale=:log10)
scatter!(monomer_concentration, BGE[2,5,:], xscale=:log10)

scatter(monomer_concentration, BGE[1,6,:], xscale=:log10)
scatter!(monomer_concentration, BGE[2,6,:], xscale=:log10)

scatter(monomer_concentration, r[1,1,:], xscale=:log10)
scatter!(monomer_concentration, r[2,1,:], xscale=:log10)

scatter(monomer_concentration, m_E[1,1,:], xscale=:log10)
scatter!(monomer_concentration, m_E[2,1,:], xscale=:log10)

scatter(monomer_concentration, m_E[1,2,:], xscale=:log10)
scatter!(monomer_concentration, m_E[2,2,:], xscale=:log10)

scatter(monomer_concentration, m_E[1,3,:], xscale=:log10)
scatter!(monomer_concentration, m_E[2,3,:], xscale=:log10)

scatter(monomer_concentration, m_E[1,4,:], xscale=:log10)
scatter!(monomer_concentration, m_E[2,4,:], xscale=:log10)

scatter(monomer_concentration, m_E[1,5,:], xscale=:log10)
scatter!(monomer_concentration, m_E[2,5,:], xscale=:log10)

scatter(monomer_concentration, m_E[1,6,:], xscale=:log10)
scatter!(monomer_concentration, m_E[2,6,:], xscale=:log10)

scatter(monomer_concentration, x[1,2,:], xscale=:log10)
scatter!(monomer_concentration, x[2,2,:], xscale=:log10)
