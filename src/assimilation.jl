abstract type AbstractMonomer end
abstract type AbstractGenomicData end
abstract type AbstractAssimilation end
abstract type AbstractMetabolism end
abstract type AbstractSoil end

@units @description struct BaseMonomer <: AbstractMonomer
    D_S::Array{Float64,1}     | m^2/s   | "Aqueous diffusivity"
    N_C::Array{Float64,1}               | "Number of C atoms"
    N_N::Array{Float64,1}               | "Number of N atoms"
    y_DE::Array{Float64,1}    | mol/mol | "Assimilation yield"
end

@units @description mutable struct PAssimilation{T<:AbstractGenomicData, D<:AbstractMonomer, G<:AbstractSoil, H<:AbstractMetabolism} <: AbstractAssimilation
   Isolate::T                           | "Isolate defined by BaseGenome"
   Monomer::D                           | "Monomer defined by BaseMonomer"
   Soil::G                              | "Soil properties defined by BaseSoil"
   PMetabolismC::H                      | "C metabolism"
   ρ_p::Array{Float64,2}                | "Transporter density on cell surface"
   N_SB::Array{Float64,2}  |mol/mol     | "Number of uptake sites"
   K_D_0::Array{Float64,2} |mol/m^3     | "Reference affinity constant"
   K_D::Array{Float64,2}   |mol/m^3     | "Affinity constant"
   #
   function PAssimilation{T, D, G, H}(Isolate::T, Monomer::D, Soil::G, PMetabolismC::H) where {T<:AbstractGenomicData, D<:AbstractMonomer, G<:AbstractSoil, H<:AbstractMetabolism}
      ρ_p = transporter_density(genome_size_to_cell_volume(Isolate.L_DNA), Isolate.min_gentime, PMetabolismC.k_E, PMetabolismC.α, PMetabolismC.y_EM, PMetabolismC.y_EV, Isolate.transporter_distr)
      N_SB = monomer_uptake_sites(genome_size_to_cell_volume(Isolate.L_DNA), ρ_p)
      K_D_0 = specific_reference_affinity(ρ_p, genome_size_to_cell_volume(Isolate.L_DNA), Monomer.D_S)
      K_D = affinity_constant(ρ_p, genome_size_to_cell_volume(Isolate.L_DNA), Monomer.D_S, Soil.pct_sand, Soil.pct_clay, Soil.s_sat, Soil.N_cell)
      new(Isolate, Monomer, Soil, PMetabolismC, ρ_p, N_SB, K_D_0, K_D)
   end
end

PAssimilation(Isolate::T, Monomer::D, Soil::G, PMetabolismC::H) where {T<:AbstractGenomicData, D<:AbstractMonomer, G<:AbstractSoil, H<:AbstractMetabolism} = PAssimilation{T, D, G, H}(Isolate::T, Monomer::D, Soil::G, PMetabolismC::H)


function specific_reference_affinity(ρ_p::Array{Float64,2}, V_c::Array{Float64,1}, D_S::Array{Float64,1})
    V_S = 100.0                   # s
    N_A = 6.022e23                # Avogrado constant [1/mol]
    r_p = 1e-9                    # 1 μm
    r_c = (3*V_c/(4*pi)).^(1/3)
    #
    p_inv = zeros(size(ρ_p,1), size(ρ_p,2))
    for i in 1:size(ρ_p,2)
        p_inv[:,i] = @. (4*ρ_p[:,i]*r_c[i]^2/r_p + pi*r_c[i])/(4*ρ_p[:,i]*r_c[i]^2/r_p)
    end
    #
    K_SC_0 = zeros(size(ρ_p,1), size(ρ_p,2))
    for i in 1:size(ρ_p,2)
        K_SC_0[:,i] = @. (V_S*ρ_p[:,i]*r_c[i])/(pi*D_S*r_p^2*N_A)*p_inv[:,i]
    end
    #
    any(x->x==true, isnan.(K_SC_0)) ? throw(DomainError("NaN in DEBJulia.ECA_kinetics!"))  : return K_SC_0
end

function affinity_constant(ρ_p::Array{Float64,2}, V_c::Array{Float64,1}, D_S::Array{Float64,1}, pct_sand::Float64, pct_clay::Float64, s_sat::Float64, N_cell::Array{Float64,1})
    V_S = 100.0
    N_A = 6.022e23
    ϕ_w, τ_g, τ_w, δ  = soil_affinity_properties(pct_sand, pct_clay, s_sat)
    K_S_0 = specific_reference_affinity(ρ_p, V_c, D_S)
    r_c = (3*V_c/(4*pi)).^(1/3)
    r_p = 1e-9
    k1_s = zeros(size(ρ_p,1), size(ρ_p,2))
    for i in 1:size(ρ_p,2)
        k1_s[:,i] = @. V_S*ρ_p[:,i]*4*r_c[i]^2/(K_S_0[:,i]*r_p^2)
    end
    #k1_s  = @. V_S*ρ_p*4*r_c^2/(K_S_0*r_p^2)
    r_m = @. r_c*(80.0*N_cell)^(1/3)
    v_m = @. 4/3*pi*r_m^3
    #
    kapvsi = zeros(size(ρ_p,1), size(ρ_p,2))
    for i in 1:size(ρ_p,2)
        kapvsi[:,i] = @. (δ/(D_S*r_m[i]*(r_m[i]+δ)) + 1/(D_S*τ_w*ϕ_w*(r_m[i]+δ)+1.e-20))*v_m[i]/(4.0*pi)
    end
    #
    K = zeros(size(ρ_p,1), size(ρ_p,2))
    for i in 1:size(ρ_p, 2)
        K[:,i] =  @. K_S_0[:,i]*(1 + (kapvsi[:,i])*k1_s[:,i]*N_cell[i]/(N_A*v_m[i]))
    end
    any(x->x==true, isnan.(K)) ? throw(DomainError("NaN in DEBJulia.ECA_kinetics!"))  : return K
end


function monomer_uptake_sites(V_c::Array{Float64,1}, ρ_p::Array{Float64,2})
    n_0 = cell_volume_to_cell_number_density(V_c)
    r_p = 1e-9                    # 1 μm
    r_c = (3*V_c/(4*pi)).^(1/3)
    N_SB = zeros(size(ρ_p,1), size(ρ_p,2))
    for i in 1:size(V_c,1)
        N_SB[:,i] = n_0[i]*4*r_c[i]^2/(r_p^2)*ρ_p[:,i]
    end
    N_SB[N_SB.<1e-16].=0.0
    return N_SB
end


function transporter_density(V_c::Array{Float64,1}, N_p::Array{Float64,2})
    r_p = 1e-9                    # 1 μm
    r_c = (3*V_c/(4*pi)).^(1/3)
    ρ_p = zeros(size(N_p,1), size(N_p,2))
    for i in 1:size(V_c,1)
        ρ_p[:,i] = N_p[:,i]*r_p^2/(4*r_c[i]^2)
    end
    return ρ_p
end


function transporter_density(V_c::Array{Float64,1}, min_gentime::Array{Float64,1}, k_E::Array{Float64,1}, α::Array{Float64,1}, y_EM::Array{Float64,1}, y_EV::Array{Float64,1}, genome_distr::Array{Float64,2})
    N_C = 10
    V_S = 100.0
    n_0 = cell_volume_to_cell_number_density(V_c)
    min_gentime = min_gentime*60*60   # in seconds
    gmax = log(2)./(min_gentime) # max. specific growth rate [1/s]
    k_M = cell_volume_to_specific_maintenance_rate(V_c)
    N_p = @. (gmax*(1+α)*y_EV + y_EM*k_M)/(N_C*V_S*n_0 - gmax*V_S*n_0*N_C/k_E)
    #
    closure   = zeros(size(genome_distr))
    for i in 1:size(closure,2)
        closure[:,i] = N_p[i]*genome_distr[:,i]./max(1e-8,sum(genome_distr[:,i]))
    end
    ceil_closure = ceil.(closure)
    ceil_closure[ceil_closure.==0.0].=1e-8
    #
    r_p = 1e-9                    # 1 μm
    r_c = (3*V_c/(4*pi)).^(1/3)
    ρ_p = zeros(size(ceil_closure,1), size(ceil_closure,2))
    for i in 1:size(V_c,1)
        ρ_p[:,i] = ceil_closure[:,i]*r_p^2/(4*r_c[i]^2)
    end
    return ρ_p
end
