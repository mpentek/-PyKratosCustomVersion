from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
#print(sys.path)

from numpy import *
from pyKratos import *

# add variables to be allocated from the list in variables.py
solution_step_variables = [
    DISPLACEMENT_X,     #Verschiebung
    DISPLACEMENT_Y,
    IS_LAGRANGIAN,
    EXTERNAL_FORCE_X,
    EXTERNAL_FORCE_Y,   
]



### ERSTELLT Modell
buffer_size = 1  # store current step and 2 in the past       
model_part = model_part.ModelPart(buffer_size, solution_step_variables)

zero_based_indices_for_nodes = True
GiDIO = gid_io.GidIO("gid_out",zero_based_indices_for_nodes)

GiDIO.ReadModelPart(model_part) ## einlesen

print(model_part)

## Static Scheme
scheme = static_scheme.StaticScheme(model_part)

# Builder And Solver
builder_and_solver = builder_and_solver.BuilderAndSolver(
    model_part, scheme)

## Solving Strategy
strategy = solving_strategy.SolvingStrategy(
    model_part, scheme, builder_and_solver)
strategy.Initialize()
strategy.Solve()

## Output
GiDIO.WriteMesh(model_part,"outmesh")
time = 0  
GiDIO.WriteNodalResults(DISPLACEMENT,model_part.NodeIterators(), time)

import plot_contour
plot_contour.PlotContour(model_part.NodeIterators(), DISPLACEMENT_X, "disp_DISPLACEMENT_X.png" )
plot_contour.PlotContour(model_part.NodeIterators(), DISPLACEMENT_Y, "disp_DISPLACEMENT_Y.png" )

import plot_system
scale = 100
plot_system.PlotSystem(model_part.NodeIterators(), model_part.ElementIterators(), DISPLACEMENT, "plot_DISPLACEMENT.png", scale)

# import plot_internal_forces
# plot_internal_forces.PlotInternalForces(model_part.ElementIterators())

print('finished')