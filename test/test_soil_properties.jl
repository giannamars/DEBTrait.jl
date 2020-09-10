
using DEBTrait, Test

pct_sand = 70.0
pct_clay = 30.0
s_sat    = 0.1
N_cell   = 10.0*ones(1)
Soil = BaseSoil(pct_sand, pct_clay, s_sat, N_cell)

chb, sat, psisat, ksat = DEBTrait.clapp_hornberger_par(Soil.pct_sand, Soil.pct_clay)
@test chb ≈ 8.76 atol=1e-6
@test sat ≈ 0.4008 atol=1e-6
@test psisat ≈ -9.1833 atol=1e-4
@test ksat ≈ 0.01085 atol=1e-4

psi, dpsidvsm = DEBTrait.cosby_psi(Soil.s_sat, psisat, sat, chb)
@test psi  ≈ -1e9 atol=1e-6
@test dpsidvsm ≈ 218562874251.49704 atol=1e-6

theta = s_sat*sat # phi_w
epsi = sat-theta
taug, tauw = DEBTrait.moldrup_tau(sat, epsi, theta, chb)
@test taug  ≈ 0.3479 atol=1e-4
@test tauw ≈ 0.000481 atol=1e-6

Dwpsi = DEBTrait.cosby_Dwpsi(s_sat, psisat, ksat, sat, chb)
@test Dwpsi ≈ 7.163253899338254e-19 atol=1e-6

s_sat = 0.10101010101010102
theta, τ_g, τ_w, film = DEBTrait.soil_affinity_properties(Soil.pct_sand, Soil.pct_clay, s_sat)
@test theta ≈ 0.04048 atol=1e-5
@test τ_g ≈ 0.347412 atol=1e-6
@test τ_w ≈ 0.000496 atol=1e-6
@test film ≈ 1e-7 atol=1e-6
