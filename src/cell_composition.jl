abstract type AbstractGenomicData end

@units @description @with_kw struct BaseGenome <: AbstractGenomicData
    L_DNA::Array{Float64,1} = 1e6*ones(1)                      | bp    | "Genome size"
    rRNA::Array{Float64,1} = 10.0*ones(1)                              | "rRNA copy number"
    min_gentime::Array{Float64,1} = 1.0*ones(1)                | hr    | "Minimum generation time"
    transporter_distr::Array{Float64,2} = 1.0*ones(6,1)                | "Relative transporter gene frequency in genome"
end

function genome_size_to_cell_volume(L_DNA::Array{Float64,1})
    # Kempes et al. (2016), Eq. 8
    v_N = 1.47e-27
    V_DNA = v_N*L_DNA
    D_0 = 3e-17
    β_D = 0.21
    V_c = (V_DNA/D_0).^(1/β_D) # m^3
end

function cell_volume_to_protein_volume(V_c::Array{Float64,1})
    # Kempes et al. (2016), Eq. 12
    P_0 = 3.42e-7
    β_p = 0.70
    V_p = P_0*V_c.^β_p # m^3
end

function cell_volume_to_ribosome_volume(V_c::Array{Float64,1})
    # Kempes et al. (2016), Eq. 13
    l_p = 975
    V_p = cell_volume_to_protein_volume(V_c)
    m_p = 5.81e-20
    d_p = 1.37e6
    v_p = m_p/d_p
    N_p = V_p/v_p
    ϕ = 6.20e-5
    η = 6.20e-5
    l_r = 4566
    r_r = 63
    μ_0 = 4e7
    β_B = 1.64
    μ = μ_0*V_c.^(β_B-1)
    N_r = @. l_p*N_p*(ϕ/μ + 1)/(r_r/μ - l_r*(η/μ + 1))
    v_r = 3.04e-24
    V_r = N_r*v_r # m^3
end

function cell_volume_to_mRNA_volume(V_c::Array{Float64,1})
    # Kempes et al. (2016), Eq. 18
    V_r = cell_volume_to_ribosome_volume(V_c)
    v_r = 3.04e-24
    N_r = V_r/v_r
    v_mRNA = 1.43e-24
    n_mRNA = 1.08
    V_mRNA = v_mRNA*n_mRNA*N_r # m^3
end

function cell_volume_to_tRNA_volume(V_c::Array{Float64,1})
    # Kempes et al. (2016), Eq. 17
    V_r = cell_volume_to_ribosome_volume(V_c)
    v_r = 3.04e-24
    N_r = V_r/v_r
    v_tRNA = 3.10e-26
    n_tRNA = 9.3
    V_tRNA = v_tRNA*n_tRNA*N_r # m^3
end

function cell_volume_to_envelope_volume(V_c::Array{Float64,1})
    # Kempes et al. (2016), Eq. 19
    r_c = (3*V_c/(4*π)).^(1/3)
    r_env = 40.5e-9      # 27.5 if gram +
    p_p = 0.149
    V_env = @. (V_c - 4/3*π*(r_c - r_env).^3)*(1-p_p) # m^3
end

function genome_size_to_rRNA_copy_number(L_DNA::Array{Float64,1})
    # Roller et al. (2016) - Supp. Fig. 1
    rRNA  = @. ceil(2^(L_DNA./1e6 - 2.8)/0.66)
end

@units @description struct PCellComposition{T<:AbstractGenomicData}
     Isolate::T                   | "Isolate defined by BaseGenome"
     V_cell::Array{Float64,1}  | m^3 | "Total cell volume"
     V_p::Array{Float64,1}     | m^3 | "Protein volume"
     V_r::Array{Float64,1}     | m^3 | "Ribosome volume"
     V_mRNA::Array{Float64,1}  | m^3 | "mRNA volume"
     V_tRNA::Array{Float64,1}  | m^3 | "tRNA volume"
     V_env::Array{Float64,1}   | m^3 | "cell envelope volume"

     function PCellComposition{T}(Isolate::T) where T<:AbstractGenomicData
         V_cell = genome_size_to_cell_volume(Isolate.L_DNA)
         V_p    = cell_volume_to_protein_volume(V_cell)
         V_r    = cell_volume_to_ribosome_volume(V_cell)
         V_mRNA = cell_volume_to_mRNA_volume(V_cell)
         V_tRNA = cell_volume_to_tRNA_volume(V_cell)
         V_env  = cell_volume_to_envelope_volume(V_cell)
         new(Isolate, V_cell, V_p, V_r, V_mRNA, V_tRNA, V_env)
     end
end


function cell_volume_to_dry_mass(V_c::Array{Float64,1})
    v_0 = 3.0e-17
    β_D = 0.21
    V_DNA = v_0*V_c.^β_D
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c)
    V_mRNA = cell_volume_to_mRNA_volume(V_c)
    V_tRNA = cell_volume_to_tRNA_volume(V_c)
    V_env = cell_volume_to_envelope_volume(V_c)
    d_DNA = 2e6
    d_p = 1.37e6
    d_r = 1.79e6
    d_mem = 1.05e6
    d_RNA = 2e6
    m_d = d_DNA*V_DNA + d_p*V_p + d_r*V_r + d_mem*V_env + d_RNA*V_mRNA + d_RNA*V_tRNA  # g
end

function cell_volume_to_cell_number_density(V_c::Array{Float64,1})
    # assuming 47% of dry mass is C (BNID 100649)
    m_d = cell_volume_to_dry_mass(V_c)
    M_C = 12.0111
    N_A = 6.022e23
    λ_B = M_C./(0.47*m_d*N_A)    # mol cells/mol C
end

function cell_volume_to_dry_mass_baseline(V_c::Array{Float64,1})
    v_0 = 3.0e-17
    β_D = 0.21
    V_DNA = v_0*V_c.^β_D
    V_env = cell_volume_to_envelope_volume(V_c)
    d_DNA = 2e6
    d_mem = 1.05e6
    m_d = d_DNA*V_DNA + d_mem*V_env  # g
end


function cell_volume_to_dry_mass_growth(V_c::Array{Float64,1})
    V_r = cell_volume_to_ribosome_volume(V_c)
    V_mRNA = cell_volume_to_mRNA_volume(V_c)
    V_tRNA = cell_volume_to_tRNA_volume(V_c)
    d_r = 1.79e6
    d_RNA = 2e6
    m_d =  d_r*V_r + d_RNA*V_mRNA + d_RNA*V_tRNA  # g
end


function reserve_yield_constant(p::PCellComposition)
    m_cell = cell_volume_to_dry_mass(p.V_cell)
    m_baseline = cell_volume_to_dry_mass_baseline(p.V_cell)
    z = m_baseline./m_cell
    X_z = [0.43, 0.143]
    m_growth = cell_volume_to_dry_mass_growth(p.V_cell)
    g = m_growth./m_cell
    X_g = [0.37, 0.153]
    y_EV = zeros(size(X_z,1), size(z,1))
    for i in 1:size(z,1)
            y_EV[:,i] = @. g[i]*X_g/(z[i]*X_z)
    end
    return y_EV
end
