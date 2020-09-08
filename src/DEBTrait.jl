module DEBTrait

using Parameters, FieldMetadata
import FieldMetadata: @default, @description, @units, @bounds, @logscaled, @flattenable, @plottable, @selectable,
                      default, description, units, bounds, logscaled, flattenable, plottable, selectable
using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ, R
using Roots

include("metabolism.jl")
export PMetabolism, growth!, growth_production

# Write your package code here.

end
