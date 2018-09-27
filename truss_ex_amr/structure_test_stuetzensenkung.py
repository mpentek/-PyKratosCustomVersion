from __future__ import print_function, absolute_import, division 
import sys
sys.path.append("..")
#print(sys.path)

from numpy import *
from pyKratos import *

# add variables to be allocated from the list in variables.py
solution_step_variables = [
    DISPLACEMENT_X,    
    DISPLACEMENT_Y,
    IS_LAGRANGIAN,
    EXTERNAL_FORCE_X,
    EXTERNAL_FORCE_Y,   
]

#      3 - - 4
#     / \   /  \
#   0 - - 1 - - 2
#         |
#         5

property_list = { 
    0: {YOUNG_MODULUS: 210e9,   # in N/m^2
        SECTION_TYPE:  0.001885, # in m^2
        }
}

node_list = {
    0: array([0.0, 0.0]), # Koordinaten in m
    1: array([5, 0]),
    2: array([10, 0]),
    3: array([2.5, -4.33]),
    4: array([7.5, -4.33]),
    5: array([5, 3]),
}

element_connectivities = {
    1: [0, [0, 1]],     ## 0 is property pointer
    2: [0, [1, 2]],     
    3: [0, [0, 3]],
    4: [0, [3, 1]],
    5: [0, [1, 4]],
    6: [0, [4, 2]],
    7: [0, [3, 4]],
    8: [0, [1, 5]],
}

condition_connectivities = {
    1: [0, [1]],
}

nodal_values = {        #Randbedingungen
DISPLACEMENT_X: [
        [0, True, 0.0], #[nodeID, fixity, imposed value in m]
        [2, True, 0.0],
        [5, True, 0.0],
    ],
DISPLACEMENT_Y: [
    [0, True, 0.0],
    [2, True, 0.0],
    [5, True, 0.2],
    ],
EXTERNAL_FORCE_Y: [
    [1, True, 000],  #[nodeID, fixity, imposed value in N]
],
}

# Model Part
buffer_size = 1     
model_part = model_part.ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements("truss_element_linear", element_connectivities)
model_part.AddNodalConditions("point_condition_2d", condition_connectivities)
model_part.AddNodalValues(nodal_values)

print(model_part)

# Static Scheme
scheme = static_scheme.StaticScheme(model_part)

# Builder And Solver
builder_and_solver = builder_and_solver.BuilderAndSolver(model_part, scheme)

# Solving Strategy
strategy = solving_strategy.SolvingStrategy(model_part, scheme, builder_and_solver)
strategy.Initialize()

strategy.Solve()


# Output
zero_based_indices_for_nodes = True
GiDIO = gid_io.GidIO("gid_out",zero_based_indices_for_nodes)
GiDIO.WriteMesh(model_part,"outmesh")

time = 0  
GiDIO.WriteNodalResults(DISPLACEMENT,model_part.NodeIterators(), time)

# import plot_contour
# plot_contour.PlotContour(model_part.NodeIterators(), DISPLACEMENT_X, "disp_DISPLACEMENT_X.png" )
# plot_contour.PlotContour(model_part.NodeIterators(), DISPLACEMENT_Y, "disp_DISPLACEMENT_Y.png" )

import plot_system
scale = 10
plot_system.PlotSystem(model_part.NodeIterators(), model_part.ElementIterators(), DISPLACEMENT, "plot_DISPLACEMENT.png", scale)

# import plot_internal_forces
# plot_internal_forces.PlotInternalForces(model_part.ElementIterators())

print('finished')