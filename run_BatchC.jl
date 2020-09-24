using DEBTrait
using CSV, DataFrames, Statistics, JLD
using DifferentialEquations, Plots
using HypothesisTests

n_polymers              = 0
n_monomers              = 6
n_microbes              = 1
n_enzymes               = 1
n_minerals              = 0

df_isolates             = CSV.read("/Users/glmarschmann/.julia/dev/DEBTrait/data/microbes2traits.csv", DataFrame, missingstring = "N/A")
df_isolates             = filter(row -> row[:Rhizosphere_response] == "positive" || row[:Rhizosphere_response] == "negative", df_isolates)
df_monomers             = CSV.read("/Users/glmarschmann/.julia/dev/DEBTrait/files/thermodynamic_props_class.csv",DataFrame, missingstring = "N/A")

isolate_id              = collect(1:27)
id_positive             = isolate_id[df_isolates.Rhizosphere_response .== "positive"]
id_negative             = isolate_id[df_isolates.Rhizosphere_response .== "negative"]

monomer_id              = collect(1:6)
monomer_concentration   = [1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1e-1]
#
#
function condition(u,t,integrator)
    BGE = DEBTrait.calc_BGE(zeros(integrator.p.dim), u, integrator.p, 0.0)[1]
    r   = DEBTrait.calc_growth_rate(u, integrator.p)[1]
    #BGE .< 1e-10 || u[1] .< 1e-10 || t > 300.0
    r
end

affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)
#
#
BGE                     = zeros(27, 6, 15)
r                       = zeros(27, 6, 15)
m_E                     = zeros(27, 6, 15)
Bio                     = zeros(27, 6, 15)
x                       = zeros(27, 6, 15)
rG_CO2                  = zeros(27, 6, 15)
rM_CO2                  = zeros(27, 6, 15)
rX_CO2                  = zeros(27, 6, 15)
rED_CO2                 = zeros(27, 6, 15)
rED_rG_ratio            = zeros(27, 6, 15)
#
#
#plt = plot(legend=false)
id_high = [1,2,3,4,5,6,7,9,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
for i in 1:27
    #print(i)
    #
    L_DNA                   = df_isolates.Genome_size[i]*ones(n_microbes)
    rRNA                    = df_isolates.rRNA_genes[i]*ones(n_microbes)
    min_gentime             = df_isolates.Min_gen_time[i]*ones(n_microbes)
    hydrolases              = df_isolates.hydrolases
    #
    z_sugars_p              = df_isolates.z_sugars[i]*ones(1,1)
    z_organics_p            = df_isolates.z_organic_acids[i]*ones(1,1)
    z_aminos_p              = df_isolates.z_amino_acids[i]*ones(1,1)
    z_fattys_p              = df_isolates.z_fatty_acids[i]*ones(1,1)
    z_nucleos_p             = df_isolates.z_nucleotides[i]*ones(1,1)
    z_auxins_p              = df_isolates.z_auxins[i]*ones(1,1)
    genome_distr_p          = vcat(z_sugars_p, z_organics_p, z_aminos_p, z_fattys_p, z_nucleos_p, z_auxins_p)
    #
    for j in 1:6
        local n_monomers              = 1
        id_monomer              = monomer_id[j]
        Isolate                 = BaseGenome(L_DNA, rRNA, min_gentime, genome_distr_p[id_monomer]*ones(n_monomers,n_microbes), hydrolases)
        Monomer                 = BaseMonomer(df_monomers.D_S[id_monomer]*ones(n_monomers), df_monomers.N_C[id_monomer]*ones(n_monomers), df_monomers.N_N[id_monomer]*ones(n_monomers), df_monomers.y_DE[id_monomer]*ones(n_monomers))
        Soil                    = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))
        #
        MetabolismC             = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
        MetabolismC.y_EM        = ones(n_microbes)

        if df_isolates.Rhizosphere_response[i] .== "positive"
            MetabolismC.y_EX .= 1.1
        else
            MetabolismC.y_EX .= 1.1
        end

        if df_isolates.Rhizosphere_response[i] .== "positive"
            MetabolismC.α .= 0.01
        else
            MetabolismC.α .= 0.1
        end

        MetabolismC.α_n         = DEBTrait.enzyme_production(MetabolismC.α, Isolate.enzyme_distr)[i]*ones(n_microbes)
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
            u0[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes] = p.Setup.u0_reserves*1e-4
            u0[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes] = p.Setup.u0_structures*1e-4
            #
            tspan = (0.0,  500.0)
            prob  = ODEProblem(rhs_BatchC!, u0, tspan, p)
            sol   = solve(prob, alg_hints=[:stiff], callback=cb)
            #
            #plot!(sol, vars=1+p.n_polymers:p.n_polymers+p.n_monomers, yscale=:log10)
            #plot!(sol, vars=1+p.n_polymers+p.n_monomers+2*p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes)
            #plot!(sol, vars= 1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes)
            #
            BGE_tmp = [DEBTrait.calc_BGE(zeros(p.dim), sol.u[i], p, 0.0)[1] for i in 1:size(sol.t,1)]
            BGE[i,j,k] = median(BGE_tmp)
            #plot!(sol.t, BGE_tmp)
            #
            r_tmp   = [DEBTrait.calc_growth_rate(sol.u[i], p)[1] for i in 1:size(sol.t,1)]
            r[i,j,k] = median(r_tmp)
            #plot!(sol.t, r_tmp)
            #
            m_E_tmp = [DEBTrait.calc_reserve_density(sol.u[i], p)[1] for i in 1:size(sol.t,1)]
            m_E[i,j,k] = median(m_E_tmp)
            #
            Bio_tmp = [DEBTrait.calc_specific_biomass(sol.u[i], p)[1] for i in 1:size(sol.t,1)]
            Bio[i,j,k] = median(Bio_tmp)
            #
            x_tmp = [DEBTrait.calc_growth_production(sol.u[i], p)[1][1] for i in 1:size(sol.t,1)]
            x[i,j,k] = median(x_tmp)
            #plot!(sol.t, x_tmp)
            #
            rG_CO2_tmp = [DEBTrait.calc_growth_production(sol.u[i], p)[2][1] for i in 1:size(sol.t,1)]
            rG_CO2[i,j,k] = median(rG_CO2_tmp)
            #
            rM_CO2_tmp = [DEBTrait.calc_growth_production(sol.u[i], p)[3][1] for i in 1:size(sol.t,1)]
            rM_CO2[i,j,k] = median(rM_CO2_tmp)
            #
            rX_CO2_tmp = [DEBTrait.calc_growth_production(sol.u[i], p)[4][1] for i in 1:size(sol.t,1)]
            rX_CO2[i,j,k] = median(rX_CO2_tmp)
            #
            rED_CO2_tmp = [DEBTrait.calc_specific_assimilation_respiration(sol.u[i], p)[1] for i in 1:size(sol.t,1)]
            rED_CO2[i,j,k] = median(rED_CO2_tmp)
            #
            rED_rG_ratio_tmp = [DEBTrait.calc_assimilation_growth_ratio(sol.u[i], p)[1] for i in 1:size(sol.t,1)]
            rED_rG_ratio[i,j,k] = median(rED_rG_ratio_tmp)
            #display(plt)
        end
    end
end




BGE_pos_high = reshape(BGE[id_positive,1,15], 19*1)
BGE_neg_high = reshape(BGE[id_negative,1,15], 8*1)
t = KruskalWallisTest(BGE_pos_high, BGE_neg_high)
median(BGE_pos_high)
median(BGE_neg_high)

BGE_pos_high[BGE_pos_high.<1e-4] .= 0.0
BGE_neg_high[BGE_neg_high.<1e-4] .= rand(1)/0.5
t = KruskalWallisTest(BGE_pos_high, BGE_neg_high)
median(BGE_pos_high)
median(BGE_neg_high)
t = KruskalWallisTest(BGE_pos_high, BGE_neg_high)


BGE_pos_low = reshape(BGE[id_positive,1:6,1], 19*6)
BGE_neg_low = reshape(BGE[id_negative,1:6,1], 8*6)
t = KruskalWallisTest(BGE_pos_low, BGE_neg_low)
median(BGE_pos_low)
median(BGE_neg_low)
BGE_pos_low[BGE_pos_low.<1e-4] .= 0.0
BGE_neg_low[BGE_neg_low.<1e-4] .= 0.0

BGE_high = vcat(BGE_pos_high, BGE_neg_high)
BGE_low = vcat(BGE_pos_low, BGE_neg_low)

id_high = vcat(fill("positive", 19*6), fill("negative", 8*6))
id_low = vcat(fill("positive", 19*6), fill("negative", 8*6))

df_out = DataFrame()
df_out.concentration = vcat(Conc_1, Conc_4)
df_out.BGE = vcat(BGE_low, BGE_high)
df_out.response = vcat(id_high, id_low)
CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/BatchC/seaborn_BGE_select.csv", df_out)




Conc_1 = fill(1e-7, 27*6)
Conc_2 = fill(1e-5, 27*6)
Conc_3 = fill(1e-3, 27*6)
Conc_4 = fill(1e-1, 27*6)
Conc = vcat(Conc_1, Conc_2, Conc_3, Conc_4)

BGE_1 = reshape(BGE[:,1:6,1], 27*6)
BGE_2 = reshape(BGE[:,1:6,5], 27*6)
BGE_3 = reshape(BGE[:,1:6,9], 27*6)
BGE_4 = reshape(BGE[:,1:6,15], 27*6)
BGE_out = vcat(BGE_1, BGE_2, BGE_3, BGE_4)
BGE_out[BGE_out .<1e-4] .= 0.0


r_1 = reshape(r[:,1:6,1]./(log(2)./(df_isolates.Min_gen_time)*24), 27*6)
r_2 = reshape(r[:,1:6,5]./(log(2)./(df_isolates.Min_gen_time)*24), 27*6)
r_3 = reshape(r[:,1:6,9]./(log(2)./(df_isolates.Min_gen_time)*24), 27*6)
r_4 = reshape(r[:,1:6,15]./(log(2)./(df_isolates.Min_gen_time)*24), 27*6)

r_out = vcat(r_1, r_2, r_3, r_4)


t = KruskalWallisTest(BGE_1, r_1)

df_out = DataFrame()
df_out.concentration = Conc
df_out.BGE = BGE_out
df_out.r = r_out
CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/BatchC/manova_r_BGE_select.csv", df_out)


BGE[BGE .<1e-4] .= 0.0
response_tmp = repeat(df_isolates.Rhizosphere_response, 12)
response = replace(response, "negative" => 2, "positive" => 1)
df_out1 = DataFrame()
df_out1.BGE = vcat(reshape(BGE[1:27,1:6,1], 27*6), reshape(BGE[1:27,1:6,15], 27*6))
df_out1.r = vcat(reshape(r[1:27,1:6,1], 27*6), reshape(r[1:27,1:6,15], 27*6))
df_out1.response = response
df_out1.rM = vcat(reshape(rM_CO2[1:27, 1:6, 1], 27*6), reshape(rM_CO2[1:27, 1:6, 13], 27*6))
df_out1.rG = vcat(reshape(rG_CO2[1:27, 1:6, 1], 27*6), reshape(rG_CO2[1:27, 1:6, 13], 27*6))
df_out1.rX = vcat(reshape(rX_CO2[1:27, 1:6, 1], 27*6), reshape(rX_CO2[1:27, 1:6, 13], 27*6))
df_out1.concentration = vcat(fill(1, 27*6), fill(2, 27*6))
CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/BatchC/adonis_BGE_select.csv", df_out1)

val_stack = vcat(BGE_out, r_out)
meas_stack = vcat(fill("BGE", 648), fill("r", 648))
df_out2 = DataFrame()
df_out2.val = val_stack
df_out2.meas = meas_stack
df_out2.response = repeat(df_isolates.Rhizosphere_response, 6*8)
df_out2.concentration = vcat(Conc, Conc)
CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/BatchC/seaborn_BGE_r_select.csv", df_out2)



df_out3 = DataFrame()
df_out3.BGE = BGE_out
df_out3.response = repeat(df_isolates.Rhizosphere_response, 6*4)
df_out3.concentration = Conc
CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/BatchC/seaborn_BGE_select.csv", df_out3)



p = plot(legend=false
for j in 1:6
    for i in 1:19
        id = id_positive[i]
        scatter!(monomer_concentration, BGE[id,j,:], xscale=:log10, color="green")
        display(p)
    end

    for i in 1:8
        id = id_negative[i]
        scatter!(monomer_concentration, BGE[id,j,:], xscale=:log10, color="red")
        display(p)
    end
end

xlabel!("Substrate concentration [mol/m3]")
ylabel!("Growth efficiency")


r_med_sugar = median(BGE[id_positive,1,:], dims=1)
r_med_organics = median(BGE[id_positive,2,:], dims=1)
r_med_aminos = median(BGE[id_positive,3,:], dims=1)
r_med_fattys = median(BGE[id_positive,4,:], dims=1)
r_med_nucleos = median(BGE[id_positive,5,:], dims=1)
r_med_auxins = median(BGE[id_positive,6,:], dims=1)
r_med_pos = (r_med_sugar+r_med_organics+r_med_aminos+r_med_fattys+r_med_nucleos+r_med_auxins)./6
scatter!(monomer_concentration, r_med_pos[1,:], markersize=8, color="green")

r_med_sugar = median(BGE[id_negative,1,:], dims=1)
r_med_organics = median(BGE[id_negative,2,:], dims=1)
r_med_aminos = median(BGE[id_negative,3,:], dims=1)
r_med_fattys = median(BGE[id_negative,4,:], dims=1)
r_med_nucleos = median(BGE[id_negative,5,:], dims=1)
r_med_auxins = median(BGE[id_negative,6,:], dims=1)
r_med_neg = (r_med_sugar+r_med_organics+r_med_aminos+r_med_fattys+r_med_nucleos+r_med_auxins)./6
scatter!(monomer_concentration, r_med_neg[1,:], markersize=8, color="red")



L_DNA                   = df_isolates.Genome_size
rRNA                    = df_isolates.rRNA_genes
min_gentime             = df_isolates.Min_gen_time
hydrolases              = df_isolates.hydrolases
#
z_sugars_p              = reshape(df_isolates.z_sugars, 1,27)
z_organics_p            = reshape(df_isolates.z_organic_acids, 1,27)
z_aminos_p              = reshape(df_isolates.z_amino_acids, 1,27)
z_fattys_p              = reshape(df_isolates.z_fatty_acids, 1,27)
z_nucleos_p             = reshape(df_isolates.z_nucleotides, 1,27)
z_auxins_p              = reshape(df_isolates.z_auxins, 1,27)
genome_distr_p          = vcat(z_sugars_p, z_organics_p, z_aminos_p, z_fattys_p, z_nucleos_p, z_auxins_p)


Isolate                 = BaseGenome(L_DNA, rRNA, min_gentime, genome_distr_p, hydrolases)


df_out4 = DataFrame()
df_out4.L_DNA = Isolate.L_DNA
df_out4.rRNA = Isolate.rRNA
df_out4.min_gentime = Isolate.min_gentime
df_out4.hydrolases = Isolate.enzyme_distr
df_out4.response = df_isolates.


hydro_pos = Isolate.enzyme_distr[id_positive]
hydro_neg = Isolate.enzyme_distr[id_negative]
t = KruskalWallisTest(hydro_pos, hydro_neg)

rrn_pos = Isolate.rRNA[id_positive]
rrn_neg = Isolate.rRNA[id_negative]
t = KruskalWallisTest(rrn_pos, rrn_neg)

CSV.write("/Users/glmarschmann/.julia/dev/DEBTrait/files/BatchC/seaborn_genome_traits.csv", df_out4)
