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
    ROTATION,
    IS_LAGRANGIAN,
    EXTERNAL_FORCE_X,
    EXTERNAL_FORCE_Y,   
    EXTERNAL_MOMENT, 
    LINE_LOAD_X,
    LINE_LOAD_Y,
    SURFACE_LOAD_X,
    SURFACE_LOAD_Y,
]


#      1
#     / \
#   0    2       

property_list = { 
    0: {YOUNG_MODULUS: 30000,   # in N/m^2
        #SECTION_TYPE:  0.000005,  # in m^2
        #MOMENT_INERTIA_AREA: 0.00005, # in m^4
        POISSON_RATIO: 0.2,
        }
}

node_list = {
    1: array([0, 2]), # Koordinaten in m
    2: array([0, 1]),
    3: array([0, 0]),
    4: array([2, 2]),
    5: array([2, 1]),
    6: array([2, 0]),
}

element_connectivities = {
    1: [0, [2, 4, 1]],     ## 0 is property pointer  
    2: [0, [2, 5, 4]],
    3: [0, [3, 5, 2]],
    4: [0, [3, 6, 5]],
}

point_conditions = {
    #0: [0, [1]],    
}

element_conditions = {  #surface load belong to 3 nodes
    1: [0, [2,4,1]],   
    2: [0, [2,5,4]], 
    3: [0, [3,5,2]],   
    4: [0, [3,6,5]],     
}
element_conditions_line = {  #line loads belong to 2 nodes
    5: [0, [1,4]],   
        
}


nodal_values = {        #boundary conditions
    DISPLACEMENT_X: [
        [1, True, 0.0], # first column is Id of node, second col if fixity, third is imposed value ggf. preforming in [m]
        [2, True, 0.0],
        [3, True, 0.0],
    ],
    DISPLACEMENT_Y: [
        [1, True, 0.0],
        [2, True, 0.0],
        [3, True, 0.0],
    ],
    ROTATION: [
        #[0, True, 0.0],
        #[2, True, 0.0],    
    ],

    EXTERNAL_FORCE_X: [
        #[1, True, 0],    # in N
    ],
    EXTERNAL_FORCE_Y: [
        #[1, True, 0],    # in N
    ],
    EXTERNAL_MOMENT: [
        #[1, True, 20000],    # in Nm
    ],
}

## Constant Line Loads
element_values = {
    LINE_LOAD_X: [
        #[1, True, 0], # in N
    ],
    LINE_LOAD_Y: [
        [1, True, -10], # in N
        #[1, True, 2000],
    ],
    SURFACE_LOAD_X: [
        #[1, True, 0], # in N
    ],
    SURFACE_LOAD_Y: [
        [1, True, -5], # in N
        [2, True, -5],
        [3, True, -5],
        [4, True, -5],
    ],
}

### ERSTELLT Modell
buffer_size = 1  # store current step and 2 in the past         
model_part = model_part.ModelPart(buffer_size, solution_step_variables)
model_part.AddNodes(node_list)
model_part.AddProperties(property_list)
model_part.AddElements("2d_triangular_shell", element_connectivities)
model_part.AddNodalConditions("point_condition_2d", point_conditions)
model_part.AddElementConditions("surface_load_condition_2d", element_conditions)
model_part.AddElementConditions("line_load_condition_2d", element_conditions_line)
model_part.AddNodalValues(nodal_values)
model_part.AddElementValues(element_values)

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

zero_based_indices_for_nodes = True
GiDIO = gid_io.GidIO("gid_out",zero_based_indices_for_nodes)
GiDIO.WriteMesh(model_part,"outmesh")

time = 0  
GiDIO.WriteNodalResults(DISPLACEMENT,model_part.NodeIterators(), time)

# import plot_system
# scale = 100
# plot_system.PlotSystem(model_part.NodeIterators(), model_part.ElementIterators(), DISPLACEMENT, "plot_DISPLACEMENT.png", scale)

# import plot_internal_forces
# plot_internal_forces.PlotInternalForces(model_part.ElementIterators())

print('finished')