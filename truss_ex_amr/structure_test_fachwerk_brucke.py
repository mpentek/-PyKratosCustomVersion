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


property_list = { 
    0: {YOUNG_MODULUS: 210e9,   # in N/m^2
        SECTION_TYPE:  0.001885, # in m^2
        }
}

node_list = {
    0: array([0.0, 0.0]), # Koordinaten in m
    1: array([1, 0]),
    2: array([2, 0]),
    3: array([3, 0]),
    4: array([4, 0]),
    5: array([5, 0]),
    6: array([6, 0]),
    7: array([7, 0]),
    8: array([8, 0]),
    9: array([9, 0]),
    10: array([10, 0]),
    11: array([11, 0]),
    12: array([12, 0]),
    13: array([13, 0]),
    14: array([14, 0]),
    15: array([15, 0]),
    16: array([16, 0]),
    17: array([17, 0]),
    18: array([18, 0]),
    19: array([19, 0]),
    20: array([20, 0]),
    21: array([0,  2]),
    22: array([1, 1.9]),
    23: array([2, 1.8]),
    24: array([3, 1.7]),
    25: array([4, 1.6]),
    26: array([5, 1.5]),
    27: array([6, 1.4]),
    28: array([7, 1.3]),
    29: array([8, 1.2]),
    30: array([9, 1.1]),
    31: array([10, 1.0]),
    32: array([11, 1.1]),
    33: array([12, 1.2]),
    34: array([13, 1.3]),
    35: array([14, 1.4]),
    36: array([15, 1.5]),
    37: array([16, 1.6]),
    38: array([17, 1.7]),
    39: array([18, 1.8]),
    40: array([19, 1.9]),
    41: array([20, 2.0]),

}

element_connectivities = {
    1: [0, [0, 1]],     ## 0 is property pointer
    2: [0, [1, 2]],     
    3: [0, [2, 3]],
    4: [0, [3, 4]],
    5: [0, [4, 5]],
    6: [0, [5, 6]],
    7: [0, [6, 7]],
    8: [0, [7, 8]],
    9: [0, [8, 9]],
    10: [0, [9, 10]],
    11: [0, [11, 12]],
    12: [0, [12, 13]],
    13: [0, [13, 14]],
    14: [0, [14, 15]],
    15: [0, [15, 16]],
    16: [0, [16, 17]],
    17: [0, [17, 18]],
    18: [0, [18, 19]],
    19: [0, [19, 20]],
    20: [0, [21, 22]],
    21: [0, [22, 23]],
    22: [0, [23, 24]],
    23: [0, [24, 25]],
    24: [0, [25, 26]],
    25: [0, [26, 27]],
    26: [0, [27, 28]],
    27: [0, [28, 29]],
    28: [0, [29, 30]],
    29: [0, [30, 31]],
    30: [0, [31, 32]],
    31: [0, [32, 33]],
    32: [0, [33, 34]],
    33: [0, [34, 35]],
    34: [0, [35, 36]],
    35: [0, [36, 37]],
    36: [0, [37, 38]],
    37: [0, [38, 39]],
    38: [0, [39, 40]],
    39: [0, [40, 41]],

    40: [0, [21, 1]],
    41: [0, [22, 2]],
    42: [0, [23, 3]],
    43: [0, [24, 4]],
    44: [0, [25, 5]],
    45: [0, [26, 6]],
    46: [0, [27, 7]],
    47: [0, [28, 8]],
    48: [0, [29, 9]],
    49: [0, [30, 10]],
    50: [0, [41, 19]],
    51: [0, [40, 18]],
    52: [0, [39, 17]],
    53: [0, [38, 16]],
    54: [0, [37, 15]],
    55: [0, [36, 14]],
    56: [0, [35, 13]],
    57: [0, [34, 12]],
    58: [0, [33, 11]],
    59: [0, [32, 10]],

    60: [0, [1, 22]],
    61: [0, [2, 23]],
    62: [0, [3, 24]],
    63: [0, [4, 25]],
    64: [0, [5, 26]],
    65: [0, [6, 27]],
    66: [0, [7, 28]],
    67: [0, [8, 29]],
    68: [0, [9, 30]],
    69: [0, [10, 31]],
    70: [0, [11, 32]],
    71: [0, [12, 33]],
    72: [0, [13, 34]],
    73: [0, [14, 35]],
    74: [0, [15, 36]],
    75: [0, [16, 37]],
    76: [0, [17, 38]],
    77: [0, [18, 39]],
    78: [0, [19, 40]],

    79: [0, [10, 11]],
}

condition_connectivities = {
    1: [0, [10]],
    2: [0, [5]],
    3: [0, [15]],
}

nodal_values = {        #Randbedingungen
DISPLACEMENT_X: [
        [0, True, 0.0], #[nodeID, fixity, imposed value in m]
        [20, True, 0.0],
        [21, True, 0.0],
        [41, True, 0.0],
    ],
DISPLACEMENT_Y: [
    [0, True, 0.0],
    [20, True, 0.0],
    [21, True, 0.0],
    [41, True, 0.0],
    ],
EXTERNAL_FORCE_Y: [
    [10, True, 200000],  #[nodeID, fixity, imposed value in N]
    [5, True, 0000],
    [15, True, 0000],
],
EXTERNAL_FORCE_X:   [
   
],
}

# Model Part
buffer_size = 1     
model_part = model_part.ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements("truss_element_nonlinear", element_connectivities)
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

for i in range(0,1):
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
scale = 20
plot_system.PlotSystem(model_part.NodeIterators(), model_part.ElementIterators(), DISPLACEMENT, "plot_DISPLACEMENT.png", scale)

# import plot_internal_forces
# plot_internal_forces.PlotInternalForces(model_part.ElementIterators())

print('finished')