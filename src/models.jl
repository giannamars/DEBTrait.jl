abstract type AbstractParameterSpace end

mutable struct SetupBatchC{T<:AbstractGenomicData}
    Isolate::T
    n_polymers::Int64
    n_monomers::Int64
    n_microbes::Int64
    n_enzymes::Int64
    n_minerals::Int64
    u0_monomers::Array{Float64,1}
    u0_reserves::Array{Float64,1}
    u0_structures::Array{Float64,1}
    function SetupBatchC{T}(Isolate::T, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals) where {T<:AbstractGenomicData}
        u0_monomers = zeros(n_monomers)
        u0_reserves = 0.5*1e11.*0.47*cell_volume_to_dry_mass(genome_size_to_cell_volume(Isolate.L_DNA))./12.011
        u0_structures = 0.5*1e11.*0.47*cell_volume_to_dry_mass(genome_size_to_cell_volume(Isolate.L_DNA))./12.011
        new(Isolate, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals, u0_monomers, u0_reserves, u0_structures)
    end
end

SetupBatchC(Isolate::T, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals) where {T<:AbstractGenomicData} = SetupBatchC{T}(Isolate::T, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)


@units @description mutable struct ParamsBatchC{B<:SetupBatchC, D<:AbstractMonomer, H<:AbstractMetabolism, A<:AbstractAssimilation} <:AbstractParameterSpace
    Setup::B                             | "System dimensions and initial conditions"
    Monomer::D                           | "Monomer defined by BaseMonomer"
    PMetabolismC::H                      | "C metabolism parameters"
    PAssimilation::A                     | "Assimilation parameters"
    #
    n_polymers::Int64
    n_monomers::Int64
    n_microbes::Int64
    n_enzymes::Int64
    n_minerals::Int64
    dim::Int64
    #
    y_DE::Array{Float64,1}   | mol/mol   | "Yield on assimilation"
    N_C::Array{Float64,1}                | "Number of C atoms"
    V_S::Array{Float64,2}    |mol/(mol*d)| "Transporter specific rate"
    N_SB::Array{Float64,2}   |mol/mol    | "Number of uptake sites"
    K_D:: Array{Float64,2}   |mol/m^3    | "Affinity constant"
    #
    γ_X::Array{Float64,1}    | d^-1      | "enzyme turnover rate"
    γ_V::Array{Float64,1}    | d^-1      | "microbial turnover rate"
    β::Array{Float64,1}                  | "density dependent turnover"
    f_E::Array{Float64,1}                | "fractionation of dead reserve into monomers"
    #
    I_D::Array{Float64,1}    |mol/(m^3*d)| "Monomer input"
    #
    function ParamsBatchC{B, D, H, A}(Setup::B, Monomer::D, PMetabolismC::H, PAssimilation::A) where {B<:SetupBatchC, D<:AbstractMonomer, H<:AbstractMetabolism, A<:AbstractAssimilation}
        n_polymers = Setup.n_polymers
        n_monomers = Setup.n_monomers
        n_microbes = Setup.n_microbes
        n_enzymes = Setup.n_enzymes
        n_minerals = Setup.n_minerals
        dim = n_polymers + n_monomers + 3*n_microbes + n_enzymes
        y_DE = Monomer.y_DE
        N_C  = Monomer.N_C
        V_S  = 100.0*24*60*60*ones(n_monomers,n_microbes)
        N_SB = PAssimilation.N_SB
        K_D = PAssimilation.K_D
        γ_X = 1e-4*ones(n_enzymes)
        γ_V = 1e-3*ones(n_microbes)
        β   = ones(n_microbes)
        f_E = ones(n_monomers).*1/n_monomers
        I_D = zeros(n_monomers)
        new(Setup, Monomer, PMetabolismC, PAssimilation, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals, dim, y_DE, N_C, V_S, N_SB, K_D, γ_X, γ_V, β, f_E, I_D)
    end
end

ParamsBatchC(Setup::B, Monomer::D, PMetabolismC::H, PAssimilation::A) where {B<:SetupBatchC, D<:AbstractMonomer, H<:AbstractMetabolism, A<:AbstractAssimilation} = ParamsBatchC{B, D, H, A}(Setup::B, Monomer::D, PMetabolismC::H, PAssimilation::A)

function rhs_BatchC!(du, u, p::AbstractParameterSpace, t)
    D                            = u[1+p.n_polymers:p.n_polymers+p.n_monomers]
    E                            = u[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    V                            = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    X                            = u[1+p.n_polymers+p.n_monomers+2*p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes]
    CO2                          = u[1+p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes+p.n_microbes]
    #
    r                            = growth!(zeros(p.n_microbes), E, V, p.PMetabolismC)
    x, rG_CO2, rM_CO2, rX_CO2    = growth_production(r, E, V, p.PMetabolismC)
    ECA_D                        = DEBTrait.ECA_kinetics!(zeros(p.n_monomers,p.n_microbes), D, V, p.K_D, p.V_S, p.N_SB)
    #
    F_D                          = vcat(sum(ECA_D[:, 1:p.n_microbes], dims=2)...).*p.N_C
    F_DE                         = vcat(sum((1.0./p.y_DE .- 1.0).*p.N_C.*ECA_D[:, 1:p.n_microbes], dims=1)...)
    F_DE_CO2                     = vcat(sum((1.0./p.y_DE).*p.N_C.*ECA_D[:, 1:p.n_microbes], dims=1)...)
    #
    F_X                          = sum(x.*V, dims=1)
    #
    F_ED                         = p.f_E.*sum(p.γ_V.*E.^p.β, dims=1)
    #
    @. du[1+p.n_polymers:p.n_polymers+p.n_monomers] = p.I_D - F_D + F_ED +  p.γ_X*X
    @. du[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes] = F_DE - (p.PMetabolismC.k_E - r)*E - p.γ_V*E^p.β
    @. du[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes] = r*V - p.γ_V*V^p.β
    @. du[1+p.n_polymers+p.n_monomers+2*p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes] = F_X - p.γ_X*X
    @. du[1+p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes+p.n_microbes]  = rG_CO2 + rM_CO2 + rX_CO2 + F_DE_CO2
    return du
end


function calc_BR(du, u, p, t)
    out = rhs_BatchC!(du, u, p, t)
    BR = out[1+p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes+p.n_microbes]
    E = u[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    V = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    V_cell = DEBTrait.genome_size_to_cell_volume(p.Setup.Isolate.L_DNA)
    Biomass = @. E + V
    N_cells = Biomass.*12.011*0.47./DEBTrait.cell_volume_to_dry_mass(V_cell)
    BR_norm = BR./N_cells
end


function calc_BP(du, u, p, t)
    out = rhs_BatchC!(du, u, p, t)
    BP_E = out[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    BP_V = out[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    BP = @. BP_E + BP_V
    V_cell = DEBTrait.genome_size_to_cell_volume(p.Setup.Isolate.L_DNA)
    E = u[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    V = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    Biomass = @. E + V
    N_cells = Biomass.*12.011*0.47./DEBTrait.cell_volume_to_dry_mass(V_cell)
    BP_norm = BP./N_cells
end


function calc_BGE(du,u,p,t)
    BP = calc_BP(du,u,p,t)
    BR = calc_BR(du,u,p,t)
    denom = max.(1e-15, BP.+BR)
    BGE = @. max.(1e-15, BP/denom)
end

function calc_growth_rate(u, p)
    E = u[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    V = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    r = growth!(zeros(p.n_microbes), E, V, p.PMetabolismC)
    return r
end

function calc_growth_production(u, p)
    E = u[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    V = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    r = growth!(zeros(p.n_microbes), E, V, p.PMetabolismC)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production(r, E, V, p.PMetabolismC)
    V_cell = DEBTrait.genome_size_to_cell_volume(p.Setup.Isolate.L_DNA)
    Biomass = @. E + V
    N_cells = Biomass.*12.011*0.47./DEBTrait.cell_volume_to_dry_mass(V_cell)
    return x./N_cells, rG_CO2./N_cells, rM_CO2./N_cells, rX_CO2./N_cells
end


function calc_specific_assimilation_rate(u, p)
    D = u[1+p.n_polymers:p.n_polymers+p.n_monomers]
    V = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    ECA_D = DEBTrait.ECA_kinetics!(zeros(p.n_monomers,p.n_microbes), D, V, p.K_D, p.V_S, p.N_SB)
    F_D   = vcat(sum(ECA_D[:, 1:p.n_microbes], dims=2)...).*p.N_C
    F_DE  = vcat(sum((1.0./p.y_DE .- 1.0).*p.N_C.*ECA_D[:, 1:p.n_microbes], dims=1)...)
    E = u[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    V_cell = DEBTrait.genome_size_to_cell_volume(p.Setup.Isolate.L_DNA)
    Biomass = @. E + V
    N_cells = Biomass.*12.011*0.47./DEBTrait.cell_volume_to_dry_mass(V_cell)
    return F_DE./N_cells
end


function calc_assimilation_rate(u, p)
    D = u[1+p.n_polymers:p.n_polymers+p.n_monomers]
    V = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    ECA_D = DEBTrait.ECA_kinetics!(zeros(p.n_monomers,p.n_microbes), D, V, p.K_D, p.V_S, p.N_SB)
    F_D   = vcat(sum(ECA_D[:, 1:p.n_microbes], dims=2)...).*p.N_C
    F_DE  = vcat(sum((1.0./p.y_DE .- 1.0).*p.N_C.*ECA_D[:, 1:p.n_microbes], dims=1)...)
end


function calc_specific_assimilation_respiration(u, p)
    D = u[1+p.n_polymers:p.n_polymers+p.n_monomers]
    V = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    ECA_D = DEBTrait.ECA_kinetics!(zeros(p.n_monomers,p.n_microbes), D, V, p.K_D, p.V_S, p.N_SB)
    E = u[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    V_cell = DEBTrait.genome_size_to_cell_volume(p.Setup.Isolate.L_DNA)
    Biomass = @. E + V
    N_cells = Biomass.*12.011*0.47./DEBTrait.cell_volume_to_dry_mass(V_cell)
    F_DE_CO2  = vcat(sum((1.0./p.y_DE).*p.N_C.*ECA_D[:, 1:p.n_microbes], dims=1)...)
    return F_DE_CO2./N_cells
end

function calc_assimilation_growth_ratio(u,p)
    F_DE_CO2 = calc_assimilation_rate(u,p)
    r        = calc_growth_rate(u,p)
    return F_DE_CO2./r
end


function calc_reserve_density(u, p)
    E = u[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    V = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    m_E = @. E/V
end

function calc_specific_biomass(u, p)
    E = u[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes]
    V = u[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes]
    Bio = @. E+V
    V_cell = DEBTrait.genome_size_to_cell_volume(p.Setup.Isolate.L_DNA)
    n_0 = DEBTrait.cell_volume_to_cell_number_density(V_cell)
    return Bio.*n_0
end
