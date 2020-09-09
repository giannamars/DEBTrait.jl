using DEBTrait, Test

Isolate          = BaseGenome(1e6*ones(1), 10*ones(1))
Cell_Composition = PCellComposition{BaseGenome}(Isolate)
y_EV             = DEBTrait.reserve_yield_constant(Cell_Composition)
