using DEBTrait, Test

################################################################################
# test BatchC
###################
n_polymers      = 0
n_monomers      = 1
n_microbes      = 1
n_enzymes       = 1
n_minerals      = 0
#
Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes), 1.0*ones(n_microbes), 1.0*ones(n_monomers,n_microbes), 1.0*ones(n_microbes))
Setup           = SetupBatchC(Isolate, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
Monomer         = BaseMonomer(1e-10*ones(n_monomers), 6*ones(n_monomers), 3*ones(n_monomers), 0.3*ones(n_monomers))
Soil            = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))
MetabolismC     = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
AssimilationC   = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
#
p = ParamsBatchC(Setup, Monomer, MetabolismC, AssimilationC)
p.PMetabolismC.y_EM = ones(n_microbes)
p.PMetabolismC.y_EX = ones(n_microbes)
p.PMetabolismC.α_n  = ones(n_microbes)*200.0
p.I_D = 1e-3*ones(n_monomers)
#
du = zeros(p.dim)
u  = 1e-8*ones(p.dim)
t = 0.0
#
out = rhs_BatchC!(du,u,p,t)
mass_balance = sum(out) .- p.I_D
@test mass_balance[1] <= 1e-10

###################
n_polymers      = 0
n_monomers      = 2
n_microbes      = 1
n_enzymes       = 1
n_minerals      = 0
#
Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes), 1.0*ones(n_microbes), 1.0*ones(n_monomers,n_microbes))
Setup           = SetupBatchC(Isolate, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
Monomer         = BaseMonomer(1e-10*ones(n_monomers), 6*ones(n_monomers), 3*ones(n_monomers), 0.3*ones(n_monomers))
Soil            = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))
MetabolismC     = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
AssimilationC   = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
#
p = ParamsBatchC(Setup, Monomer, MetabolismC, AssimilationC)
p.PMetabolismC.y_EM = ones(n_microbes)
p.PMetabolismC.y_EX = ones(n_microbes)
p.PMetabolismC.α   = zeros(n_microbes)
p.I_D = 1e-3*ones(n_monomers)
#
du = zeros(p.dim)
u  = 1e-8*ones(p.dim)
t = 0.0
#
out = rhs_BatchC!(du,u,p,t)
mass_balance = sum(out) .- sum(p.I_D)
@test mass_balance[1] <= 1e-10

###################
n_polymers      = 0
n_monomers      = 1
n_microbes      = 2
n_enzymes       = 1
n_minerals      = 0
#
Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes), 1.0*ones(n_microbes), 1.0*ones(n_monomers,n_microbes))
Setup           = SetupBatchC(Isolate, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
Monomer         = BaseMonomer(1e-10*ones(n_monomers), 6*ones(n_monomers), 3*ones(n_monomers), 0.3*ones(n_monomers))
Soil            = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))
MetabolismC     = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
AssimilationC   = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
#
p = ParamsBatchC(Setup, Monomer, MetabolismC, AssimilationC)
p.PMetabolismC.y_EM = ones(n_microbes)
p.PMetabolismC.y_EX = ones(n_microbes)
p.PMetabolismC.α   = zeros(n_microbes)
p.I_D = 1e-3*ones(n_monomers)
#
du = zeros(p.dim)
u  = 1e-8*ones(p.dim)
t = 0.0
#
out = rhs_BatchC!(du,u,p,t)
mass_balance = sum(out) .- sum(p.I_D)
@test mass_balance[1] <= 1e-10

#################
n_polymers      = 0
n_monomers      = 2
n_microbes      = 2
n_enzymes       = 1
n_minerals      = 0
#
Isolate         = BaseGenome(1e6*ones(n_microbes), 10*ones(n_microbes), 1.0*ones(n_microbes), 1.0*ones(n_monomers,n_microbes))
Setup           = SetupBatchC(Isolate, n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
Monomer         = BaseMonomer(1e-10*ones(n_monomers), 6*ones(n_monomers), 3*ones(n_monomers), 0.3*ones(n_monomers))
Soil            = BaseSoil(0.0, 0.0, 1.0, 1*ones(n_microbes))
MetabolismC     = PMetabolismC{BaseGenome, Array{Float64,1}}(Isolate)
AssimilationC   = PAssimilation(Isolate, Monomer, Soil, MetabolismC)
#
p = ParamsBatchC(Setup, Monomer, MetabolismC, AssimilationC)
p.PMetabolismC.y_EM = ones(n_microbes)
p.PMetabolismC.y_EX = ones(n_microbes)
p.PMetabolismC.α   = zeros(n_microbes)
p.I_D = 1e-3*ones(n_monomers)
#
du = zeros(p.dim)
u  = 1e-8*ones(p.dim)
t = 0.0
#
out = rhs_BatchC!(du,u,p,t)
mass_balance = sum(out) .- sum(p.I_D)
@test mass_balance[1] <= 1e-10

################################################################################
