module DEBTrait

using Parameters, FieldMetadata
import FieldMetadata: @default, @description, @units, @bounds, @logscaled, @flattenable, @plottable, @selectable,
                      default, description, units, bounds, logscaled, flattenable, plottable, selectable
using Unitful
Unitful.register(@__MODULE__);
@unit bp "bp" Basepair 1u"m" true;
using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

using Roots
using LinearAlgebra

include("metabolism.jl")
include("cell_composition.jl")
include("assimilation.jl")
include("supeca.jl")
include("soil_properties.jl")
include("ThermoStoichWizard.jl")
include("models.jl")



export AbstractMetabolism, PMetabolismC, PMetabolismCN, growth!, growth_production
export AbstractGenomicData, BaseGenome, PCellComposition
export AbstractMonomer, AbstractAssimilation, BaseMonomer, PAssimilation
export AbstractSoil, BaseSoil
export AbstractParameterSpace, SetupBatchC, ParamsBatchC, rhs_BatchC!


# Write your package code here.

end
