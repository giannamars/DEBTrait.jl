abstract type AbstractParameterSpace end

mutable struct SetupBatchC
    n_polymers::Int64
    n_monomers::Int64
    n_microbes::Int64
    n_enzymes::Int64
    n_minerals::Int64
end

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
        f_E = ones(n_monomers).*1/n_monomers
        I_D = zeros(n_monomers)
        new(Setup, Monomer, PMetabolismC, PAssimilation, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals, dim, y_DE, N_C, V_S, N_SB, K_D, γ_X, γ_V, f_E, I_D)
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
    F_ED                         = p.f_E.*sum(p.γ_V.*E, dims=1)
    #
    @. du[1+p.n_polymers:p.n_polymers+p.n_monomers] = p.I_D - F_D + F_ED +  p.γ_X*X
    @. du[1+p.n_polymers+p.n_monomers:p.n_polymers+p.n_monomers+p.n_microbes] = F_DE - (p.PMetabolismC.k_E - r + p.γ_V)*E
    @. du[1+p.n_polymers+p.n_monomers+p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes] = (r - p.γ_V)*V
    @. du[1+p.n_polymers+p.n_monomers+2*p.n_microbes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes] = F_X - p.γ_X*X
    @. du[1+p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes:p.n_polymers+p.n_monomers+2*p.n_microbes+p.n_enzymes+p.n_microbes]  = rG_CO2 + rM_CO2 + rX_CO2 + F_DE_CO2
    return du
end
