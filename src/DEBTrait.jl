module DEBTrait

using Parameters, FieldMetadata
import FieldMetadata: @default, @description, @units, @bounds, @logscaled, @flattenable, @plottable, @selectable,
                      default, description, units, bounds, logscaled, flattenable, plottable, selectable
using Unitful
Unitful.register(@__MODULE__);
@unit bp "bp" Basepair 1u"m" true;
using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R

using Roots

include("metabolism.jl")
include("cell_composition.jl")

export AbstractMetabolism, PMetabolismC, PMetabolismCN, growth!, growth_production
export AbstractGenomicData, BaseGenome, PCellComposition


# Write your package code here.

end
