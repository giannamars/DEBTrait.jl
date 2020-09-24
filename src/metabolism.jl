abstract type AbstractGenomicData end
abstract type AbstractMetabolism end

@units @description mutable struct PMetabolismCN{T<:AbstractGenomicData, D1<:AbstractArray, D2<:AbstractArray} <: AbstractMetabolism
   Isolate::T              | "Isolate defined by BaseGenome"
   y_EM::D1      | mol/mol | "yield of maintenance on reserve"
   y_EX::D1      | mol/mol | "yield of enzyme on reserve"
   y_EV::D1      | mol/mol | "yield of structure on reserve"
   k_E::D2       | d^-1    | "specific reserve turnover rate"
   k_M::D2       | d^-1    | "specific maintenance rate"
   α::D2                   | "fraction of growth flux used for enzyme production"
   α_n::D2                 | "normalized fraction of growth flux used for enzyme production"

   function PMetabolismCN{T, D1, D2}(Isolate::T) where {T<:AbstractGenomicData, D1<:AbstractArray, D2<:AbstractArray}
         y_EV = reserve_yield_constant(genome_size_to_cell_volume(Isolate.L_DNA))
         y_EM = similar(y_EV)
         y_EX = similar(y_EV)
         k_E  = rRNA_copy_number_to_translation_power(Isolate.L_DNA, Isolate.rRNA, Isolate.min_gentime)
         k_M  = cell_volume_to_specific_maintenance_rate(genome_size_to_cell_volume(Isolate.L_DNA))
         α    = similar(k_M)
         α_n  = enzyme_production(α, Isolate.enzyme_distr)
         new(Isolate, y_EM, y_EX, y_EV, k_E, k_M, α, α_n)
     end
end

@units @description mutable struct PMetabolismC{T<:AbstractGenomicData, D1<:AbstractArray} <: AbstractMetabolism
   Isolate::T              | "Isolate defined by BaseGenome"
   y_EM::D1      | mol/mol | "yield of maintenance on reserve"
   y_EX::D1      | mol/mol | "yield of enzyme on reserve"
   y_EV::D1      | mol/mol | "yield of structure on reserve"
   k_E::D1       | d^-1    | "specific reserve turnover rate"
   k_M::D1       | d^-1    | "specific maintenance rate"
   α::D1                   | "fraction of growth flux used for enzyme production"
   α_n::D1                 | "normalized fraction of growth flux used for enzyme production"

   function PMetabolismC{T, D1}(Isolate::T) where {T<:AbstractGenomicData, D1<:AbstractArray}
         y_EV = metabolic_partitioning(genome_size_to_cell_volume(Isolate.L_DNA))
         #y_EV = reserve_yield_constant(genome_size_to_cell_volume(Isolate.L_DNA))[1,:]
         y_EM = similar(y_EV)
         y_EX = similar(y_EV)
         k_E  = rRNA_copy_number_to_translation_power(Isolate.L_DNA, Isolate.rRNA, Isolate.min_gentime)
         k_M  = cell_volume_to_specific_maintenance_rate(genome_size_to_cell_volume(Isolate.L_DNA))
         α    = similar(k_M)
         α_n  = enzyme_production(α, Isolate.enzyme_distr)
         new(Isolate, y_EM, y_EX, y_EV, k_E, k_M, α, α_n)
     end
end


function growth!(r0::Array{Float64,1}, E::Array{Float64,1}, V::Array{Float64,1}, p::PMetabolismC)
   y_VM = p.y_EM./p.y_EV
   #
   function f(r, j)
     m_E     = max(1e-8, E[j]/V[j])
     j_EC    = m_E*(p.k_E[j] - r)
     j_EM    = p.k_M[j]*p.y_EM[j]
     jEM     = max(1e-8, min(j_EC, j_EM))
     jVM     = (j_EM - jEM)*y_VM[j]/p.y_EM[j]
     j_EG    = (j_EC - jEM)/(1.0 + p.α_n[j])
     j_G     = j_EG/p.y_EV[j]
     res     = r - (j_G - jVM)
   end
   #
   f1 = [r->f(r,j) for j in 1:size(r0,1)]
   r = Roots.find_zero.(f1, r0)
 end

 function growth!(r0::Array{Float64,1}, E::Array{Float64,2}, V::Array{Float64,1}, p::PMetabolismCN)
    y_VM = p.y_EM./p.y_EV
    #
    function f(r, j)
      m_E     = @. E[:,j]/V[j]
      j_EC    = @. m_E*(p.k_E - r)
      j_EM    = @. p.k_M*p.y_EM[:,j]
      jEM     = @. min(j_EC, j_EM)
      jVM     = @. (j_EM - jEM)*y_VM[:,j]/p.y_EM[:,j]
      j_EG    = @. (j_EC - jEM)/(1.0 + p.α_n[j])
      j_EG_V  = @. max(1e-6, j_EG./p.y_EV[:,j])
      j_G     = 1.0./(sum(1.0./j_EG_V) - 1.0./sum(j_EG_V))
      res     = r - (j_G - sum(jVM))
    end
    #
    f1 = [r->f(r,j) for j in 1:size(r0,1)]
    r = Roots.find_zero.(f1, r0)
  end

function growth_production(r::Array{Float64,1}, E::Array{Float64,1}, V::Array{Float64,1}, p::PMetabolismC)
   y_VM    = p.y_EM./p.y_EV
   m_E     = @. E/V
   j_EC    = @. m_E*(p.k_E - r)
   j_EM    = @. p.k_M*p.y_EM
   jEM     = @. min(j_EC, j_EM)
   jVM     = @. (j_EM - jEM)*y_VM/p.y_EM
   j_EG    = @. (j_EC - jEM)/(1.0 + p.α_n)
   x       = @. p.α_n*j_EG/p.y_EX
   #
   rG_CO2  = @. (1-1/p.y_EV)*j_EG*V
   rM_CO2  = @. (jEM + jVM)*V
   rX_CO2  = @. (1 - 1/p.y_EX)*p.α_n*j_EG*V
   #
   return x, rG_CO2, rM_CO2, rX_CO2
end


function growth_production(r::Array{Float64,1}, E::Array{Float64,2}, V::Array{Float64,1}, p::PMetabolismCN)
  y_VM    = p.y_EM./p.y_EV
  #
  x = zeros(size(p.α))
  j_Ex   = zeros(size(p.y_EM))
  rG_CO2 = zeros(size(p.y_EM))
  rM_CO2 = zeros(size(p.y_EM))
  rX_CO2 = zeros(size(p.y_EM))
  rEG_rej = zeros(size(p.y_EM))
  rEX_rej = zeros(size(p.y_EM))
  #
  for j in 1:size(r,1)
    m_E     = @. E[:,j]/V[j]
    j_EC    = @. m_E*(p.k_E - r)
    j_EM    = @. p.k_M*p.y_EM[:,j]
    jEM     = @. min(j_EC, j_EM)
    jVM     = @. (j_EM - jEM)*y_VM[:,j]/p.y_EM[:,j]
    j_EG    = @. (j_EC - jEM)/(1.0 + p.α_n[j])
    #
    j_Ex    = @. max(1e-8, p.α_n[j]*j_EG/p.y_EX[:,j])
    x[j]    = 1.0./(sum(1.0./j_Ex) - 1.0./sum(j_Ex))
    #x[isnan.(x)] .= 0.0
    #rG_CO2[:,j]  = @. (1-1/p.y_EV[:,j])*j_EG*V[j]
    rM_CO2[:,j]  = @. (jEM + jVM)*V[j]
    #rX_CO2[:,j]  = @. (1 - 1/p.y_EX[:,j])*p.α_n[j]*j_EG*V[j]
    rEG_rej[:,j]  = max.(0.0, (p.k_E - r).* m_E - jEM - p.y_EV[:,j]*(r[j] + sum(jVM)))
    rEX_rej[:,j] = max.(0.0, p.α_n[j].*((p.k_E - r).* m_E - jEM) - p.y_EX[:,j]*x[j])
  end
  return x, rM_CO2, rEG_rej, rEX_rej
end

function cell_volume_to_specific_maintenance_rate(V_c::Array{Float64,1})
    # Lynch & Marinov (2015), Eq. 1a
    n_0 = cell_volume_to_cell_number_density(V_c)
    V_cell = V_c./1e-18         #  V_cell in μm^3
    E_M = 0.39*V_cell.^0.88     # [1e9 ATP/(cell h)]
    E_G = 1e9*E_M/26            # 26 ATP per glucose molecule: [glucose/(cell h)]
    E_C = E_G*6*24              # [C atoms/(cell/d)]
    mol_C = E_C/(6.022e23)      # [mol C/(cell d)]
    k_M = mol_C./n_0            # 1/d
end

function reserve_yield_constant(V_c::Array{Float64,1})
    m_cell = cell_volume_to_dry_mass(V_c)
    m_baseline = cell_volume_to_dry_mass_baseline(V_c)
    z = m_baseline./m_cell
    X_z = [0.43, 0.143]
    m_growth = cell_volume_to_dry_mass_growth(V_c)
    g = m_growth./m_cell
    X_g = [0.37, 0.153]
    y_EV = zeros(size(X_z,1), size(z,1))
    for i in 1:size(z,1)
            y_EV[:,i] = @. g[i]*X_g/(z[i]*X_z)
    end
    return y_EV
end


function rRNA_copy_number_to_translation_power(rRNA::Array{Float64,1})
    log2_rna = log2.(rRNA)
    y = @. 0.8 + 0.61*log2_rna
    fg = y*24 #per day
end

function rRNA_copy_number_to_translation_power(L_DNA::Array{Float64,1}, rRNA::Array{Float64,1}, min_gentime::Array{Float64,1})
    log2_rna = log2.(rRNA)
    y = @. 0.8 + 0.61*log2_rna
    fg = y #per day
    V_c = genome_size_to_cell_volume(L_DNA)
    V_r = cell_volume_to_ribosome_volume(V_c)
    V_p = cell_volume_to_protein_volume(V_c)
    α = @. V_r./V_p
    k_E = fg.*α
    μ_max = @. log(2)/(min_gentime/24)
    k_E[k_E.>μ_max] = 0.95*μ_max[k_E.>μ_max]
    return k_E
end


function metabolic_partitioning(V_c::Array{Float64,1})
    μ_0 = 4e7
    β_B = 1.64
    μ = μ_0*V_c.^(β_B-1)
    b = 2.42e-6
    γ = @. 1/(1 + b/μ)
    return @. 1.0/γ
end


function enzyme_production(α::Array{Float64,1}, genome_distr::Array{Float64,1})
  genome_distr_norm = genome_distr./sum(genome_distr)
  closure = α.*genome_distr_norm
end
