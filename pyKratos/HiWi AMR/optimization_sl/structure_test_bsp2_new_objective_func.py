from __future__ import print_function, absolute_import, division
import scipy
import sys
sys.path.append("..")
#print(sys.path)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

#   0 - - 1 - - 2 - - 3
#      \  |  \  |  \  |
#   4 - - 5 - - 6 - - 7

node_list = {
    0: array([0.0, 0.0]), # Koordinaten in m
    1: array([5, 0]),
    2: array([10, 0]),
    3: array([15, 0]),
    4: array([0, -5]),
    5: array([5, -5]),
    6: array([10, -5]),
    7: array([15, -5]),
}

property_list = {
    0: {YOUNG_MODULUS: 200e9,   # in N/m^2
        SECTION_TYPE:  0.00005, # in m^2
        }
}

element_connectivities = {
    1: [0, [0, 1]],     # 0 is property pointer
    2: [0, [1, 2]],
    3: [0, [2, 3]],
    4: [0, [0, 5]],
    5: [0, [1, 5]],
    6: [0, [1, 6]],
    7: [0, [2, 6]],
    8: [0, [2, 7]],
    9: [0, [3, 7]],
    10: [0, [4, 5]],
    11: [0, [5, 6]],
    12: [0, [6, 7]],
}

condition_connectivities = {
    1: [0, [3]],
    2: [0, [7]],
}

nodal_values = {        #Randbedingungen
DISPLACEMENT_X: [
        [0, True, 0.0], #[nodeID, fixity, imposed value in m]
        [4, True, 0.0],
    ],
DISPLACEMENT_Y: [
    [0, True, 0.0],
    [4, True, 0.0],
    ],
EXTERNAL_FORCE_Y: [
    [3, True, -10000],  #[nodeID, fixity, imposed value in N]
    [7, True, -10000],
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
GiDIO = gid_io.GidIO("inputfile.mdpa","gid_out",zero_based_indices_for_nodes)
GiDIO.WriteMesh(model_part,"outmesh")

time = 0
GiDIO.WriteNodalResults(DISPLACEMENT,model_part.NodeIterators(), time)


# initial values used in plot
x_original_1 = []
y_original_1 = []
x_original_2 = []
y_original_2 = []

for element in model_part.ElementIterators():

    x_original_1.append(element.geometry.nodes[0].coordinates[0])
    y_original_1.append(element.geometry.nodes[0].coordinates[1])
    x_original_2.append(element.geometry.nodes[1].coordinates[0])
    y_original_2.append(element.geometry.nodes[1].coordinates[1])


################################################################################################
################################################################################################

# Initialize nodes to be modified
x0 = []
x0.extend(model_part.Nodes[1].coordinates)
x0.extend(model_part.Nodes[2].coordinates)
x0.extend(model_part.Nodes[5].coordinates)
x0.extend(model_part.Nodes[6].coordinates)
#x0.extend(model_part.Nodes[3].coordinates)

# Displacement before modification
w1_start = model_part.Nodes[7].variables[0]["DISPLACEMENT_Y"]
print('w1_initial',w1_start)

# objective function
def objective_function(x0):
    
    for u in range(len(model_part.Nodes)):
        model_part.Nodes[u].variables[0]["DISPLACEMENT_Y"]=0
        model_part.Nodes[u].variables[0]["DISPLACEMENT_X"]=0

    model_part.Nodes[1].coordinates = [x0[0],x0[1]]
    model_part.Nodes[2].coordinates = [x0[2],x0[3]]
    model_part.Nodes[5].coordinates = [x0[4],x0[5]]
    model_part.Nodes[6].coordinates = [x0[6],x0[7]]
    #model_part.Nodes[3].coordinates = [x0[8],x0[9]]

    strategy.Solve()
    w1 = model_part.Nodes[7].variables[0]["DISPLACEMENT_Y"]
   
    F=((1/w1_start)*w1)**2
    
    return F


##---   ----    ----       SCIPY      ----   ----   ----   ----
##---   ----    ----    ----   ----   ----   ----   ----   ----
import numpy as np
res=scipy.optimize.minimize(objective_function,x0,options={'maxiter':50})

print('')
print(' Optimization Results Scipy')
print('')
print('Best x:',res.x)
print('func:',res.fun)
print('nfev:',res.nfev)
print('nit:',res.nit)
norm=np.linalg.norm(res.jac)

print('gradientnorm:', norm)

##---   ----    ----    ----   ----   ----   ----   ----   ----
##---   ----    ----    ----   ----   ----   ----   ----   ----


# Initialize variables
print(' ')
print('----- Steepest Descent -----')

converge=False
j=0
epsilon=0.0001
delta = 0.001
delta_vec1 = zeros(len(x0))
delta_vec2 = zeros(len(x0))
delta_vec =zeros((len(x0),len(x0)))
for i in range(len(x0)):
    delta_vec1[i]=delta
    delta_vec[i]=delta_vec1
    delta_vec1[i]=0
    delta_vec2[i]=delta
evaluated_design_variables=[]
evaluated_results=[]
evaluated_gradientnorm=[]
count_iterates =[]
gamma=[]
gamma.append(0.01)
gamma_inner=np.zeros(50)
gamma_inner[0]=0.01

### ----------------------- ###

# function evaluation
F = objective_function(x0)   
df=zeros(len(x0))
for i in range(len(x0)):
    F_delta = objective_function(np.add(x0, delta_vec[i]))
    df[i]=-((F_delta-F)/(delta)) 
df_norm = np.linalg.norm(df)

cur_x=x0
print(' ')
print('Startwert:',cur_x,'Funktionswert:',F)
print('Iterationen:')

while   j < 75 and df_norm >= (epsilon) and not converge:
    i=0
    gamma.append(0)

###### NEWTON - LIKE    
    
    phi1 = (objective_function(np.add(cur_x, (gamma[j] + delta) * df))- objective_function(np.add(cur_x, (gamma[j]) * df))) / delta
    phi2 = (-2*objective_function(np.add(cur_x, (gamma[j]) * df))   +    objective_function(np.add(cur_x, (gamma[j] + delta) * df))    +    objective_function(np.add(cur_x, (gamma[j] - delta) * df)))  / delta**2
    
    q=phi1/phi2
    gamma_inner[i+1]=(gamma_inner[i] - float( q ))
    
    while   abs((gamma_inner[i+1])-gamma_inner[i]) >= epsilon:
            i=i+1
            
            phi1 = (objective_function(np.add(cur_x, (gamma_inner[i] + delta) * df))- objective_function(np.add(cur_x, (gamma_inner[i]) * df))) / delta
            phi2 = (-2*objective_function(np.add(cur_x, (gamma_inner[i]) * df))   +    objective_function(np.add(cur_x, (gamma_inner[i] + delta) * df))    +    objective_function(np.add(cur_x, (gamma_inner[i] - delta) * df)))  / delta**2
            
            q=phi1/phi2
            
            gamma_inner[i+1]=(gamma_inner[i] - float( q ))
            
            if i>20:
                    print('Newton linesearch failed')
                    break
    
    if gamma_inner[i+1] > 0 :                
            gamma[j]=(gamma_inner[i+1])
    else:
            pass
            
    print('cur_x:',cur_x,'F:',F)
    
    prev_x = cur_x
    prev_F = F
    prev_df = df

    cur_x = np.add(cur_x, gamma[j] *df)

    if np.linalg.norm(cur_x-prev_x)  <= epsilon:

        converge = True
        print('converge = True')

    else:
        j=j+1


    F = objective_function(cur_x)

    for i in range(len(x0)):
        F_delta = objective_function(np.add(cur_x, delta_vec[i]))
        df[i]=-((F_delta-F)/(delta))
    
    df_norm = np.linalg.norm(df)

    evaluated_design_variables.extend([prev_x])
    evaluated_results.append(prev_F)
    evaluated_gradientnorm.append(df_norm)
    count_iterates.append(j)
print(' ')
print('iterates:',j,'  best x:',prev_x,'  func:',prev_F, 'gradientnorm:',df_norm, 'Displacement',model_part.Nodes[7].variables[0]["DISPLACEMENT_Y"])
print(' ')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(count_iterates,evaluated_results,'k')
plt.ylabel('   f(x)',size=12)
plt.xlabel('iterations',size=12)

fig = plt.figure()
plt.plot(count_iterates,evaluated_gradientnorm,'k-')
plt.ylabel('   ||f '"'"'(x)||',size=12)
plt.xlabel('iterations',size=12)


##################### plot Fachwerk ################################

import plot_system
scale = 10 # scale factor of displacement
plot_system.PlotSystem(model_part.NodeIterators(), model_part.ElementIterators(), DISPLACEMENT, "plot_DISPLACEMENT.png", scale, x_original_1, y_original_1, x_original_2, y_original_2)


#########################################################################
#########################################################################
