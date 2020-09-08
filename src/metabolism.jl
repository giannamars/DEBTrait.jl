@units @description @with_kw mutable struct PMetabolism{T<:AbstractArray, D<:AbstractArray}
    y_EM::T = fill(1.0,1)      | mol/mol | "yield of maintenance on reserve"
    y_EV::T = fill(1.0,1)      | mol/mol | "yield of structure on reserve"
    y_EX::T = fill(1.0,1)      | mol/mol | "yield of enzyme on reserve"
    k_E::D  = fill(0.5,1)      | d^-1    | "specific reserve turnover rate"
    k_M::D  = fill(1.0,1)      | d^-1    | "specific maintenance rate"
    α::D    = fill(0.0,1)      |           "fraction of growth flux used for enzyme production"
end

function growth!(r0::Array{Float64,1}, E::Array{Float64,1}, V::Array{Float64,1}, p)
   y_VM = p.y_EM./p.y_EV
   #
   function f(r, j)
     m_E     = max(1e-8, E[j]/V[j])
     j_EC    = m_E*(p.k_E[j] - r)
     j_EM    = p.k_M[j]*p.y_EM[j]
     jEM     = max(1e-8, min(j_EC, j_EM))
     jVM     = (j_EM - jEM)*y_VM[j]/p.y_EM[j]
     j_EG    = (j_EC - jEM)/(1.0 + p.α[j])
     j_G     = j_EG/p.y_EV[j]
     res     = r - (j_G - jVM)
   end
   #
   f1 = [r->f(r,j) for j in 1:size(r0,1)]
   r = Roots.find_zero.(f1, r0)
 end

 function growth!(r0::Array{Float64,1}, E::Array{Float64,2}, V::Array{Float64,1}, p)
    y_VM = p.y_EM./p.y_EV
    #
    function f(r, j)
      m_E     = @. E[:,j]/V[j]
      j_EC    = @. m_E*(p.k_E - r)
      j_EM    = @. p.k_M*p.y_EM[:,j]
      jEM     = @. min(j_EC, j_EM)
      jVM     = @. (j_EM - jEM)*y_VM[:,j]/p.y_EM[:,j]
      j_EG    = @. (j_EC - jEM)/(1.0 + p.α[j])
      j_EG_V  = @. max(1e-6, j_EG./p.y_EV[:,j])
      j_G     = 1.0./(sum(1.0./j_EG_V) - 1.0./sum(j_EG_V))
      res     = r - (j_G - sum(jVM))
    end
    #
    f1 = [r->f(r,j) for j in 1:size(r0,1)]
    r = Roots.find_zero.(f1, r0)
  end

function growth_production(r::Array{Float64,1}, E::Array{Float64,1}, V::Array{Float64,1}, p)
   y_VM    = p.y_EM./p.y_EV
   m_E     = @. E/V
   j_EC    = @. m_E*(p.k_E - r)
   j_EM    = @. p.k_M*p.y_EM
   jEM     = @. min(j_EC, j_EM)
   jVM     = @. (j_EM - jEM)*y_VM/p.y_EM
   j_EG    = @. (j_EC - jEM)/(1.0 + p.α)
   x       = @. p.α*j_EG/p.y_EX
   #
   rG_CO2  = @. (1-1/p.y_EV)*j_EG*V
   rM_CO2  = @. (jEM + jVM)*V
   rX_CO2  = @. (1 - 1/p.y_EX)*p.α*j_EG*V
   #
   return x, rG_CO2, rM_CO2, rX_CO2
end


function growth_production(r::Array{Float64,1}, E::Array{Float64,2}, V::Array{Float64,1}, p)
  y_VM    = p.y_EM./p.y_EV
  #
  #x = zeros(size(p.y_EM))
  j_Ex   = zeros(size(p.y_EM))
  rG_CO2 = zeros(size(p.y_EM))
  rM_CO2 = zeros(size(p.y_EM))
  rX_CO2 = zeros(size(p.y_EM))
  rE_rej = zeros(size(p.y_EM))
  #
  for j in 1:size(r,1)
    m_E     = @. E[:,j]/V[j]
    j_EC    = @. m_E*(p.k_E - r)
    j_EM    = @. p.k_M*p.y_EM[:,j]
    jEM     = @. min(j_EC, j_EM)
    jVM     = @. (j_EM - jEM)*y_VM[:,j]/p.y_EM[:,j]
    j_EG    = @. (j_EC - jEM)/(1.0 + p.α[j])

      #x[j]    = 1/(sum(1/j_Ex) - 1/sum(j_Ex))
    j_Ex[:,j]    = @. p.α[j]*j_EG/p.y_EX[:,j]
    rG_CO2[:,j]  = @. (1-1/p.y_EV[:,j])*j_EG*V[j]
    rM_CO2[:,j]  = @. (jEM + jVM)*V[j]
    rX_CO2[:,j]  = @. (1 - 1/p.y_EX[:,j])*p.α[j]*j_EG*V[j]
    rE_rej[:,j]  = max.(0.0, (p.k_E - r).* m_E - jEM - p.y_EV[:,j]*(r[j] + sum(jVM)))
  end
  return j_Ex, rG_CO2, rM_CO2, rX_CO2, rE_rej
end
