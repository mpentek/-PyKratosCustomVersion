from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create( Id, prop, list_of_nodes):
    geom = line2d.Line2D(list_of_nodes)
    return LineLoad2DCondition(Id, prop, geom)

class LineLoad2DCondition:

    def __init__(self, Id, prop, geom):
        self.Id = Id
        self.prop = prop
        self.geometry = geom
        self.element = []

    def GetDofsPerNode(self):
        return 3

                
    def CalculateLocalSystem(self, ProcessInfo, model_part):
        list_of_nodes = self.geometry.nodes

        # find element Id for given Nodes
        for ele in model_part.Elements.values():
            if ele.geometry.nodes == list_of_nodes:
                self.element = ele             
            
        try: qx = self.element.GetSolutionStepValue(LINE_LOAD_X,0)
        except: qx = 0
        try: qy = self.element.GetSolutionStepValue(LINE_LOAD_Y,0)  
        except: qy = 0 
        
        RHS = zeros(6)  
        LHS = zeros((6,6))

        order = 3 # for constant line load order = 1 is enough
        [gpc, weights] = self.geometry.GaussPoints(order)
        [Ns_bending, dderivatives, ddderivatives] = self.geometry.ShapeFunctionsBending(gpc)
        [Ns_linear, derivatives] = self.geometry.ShapeFunctionsLinear(gpc)

        number_of_gauss = len(gpc)     
        
         #loop over all gauss points                
        for gauss in range(0, number_of_gauss):
            GPW = weights[gauss]

            # local load vector from linear element
            for i,j in zip( range(0, 2), array((0, 3)) ):
                RHS[j] += qx * Ns_linear[gauss][i] * GPW #GPW = WGP * L     

            # local load vector from bending
            for i,j in zip( range(0, 4), array((1, 2, 4, 5)) ):
                RHS[j] += qy * Ns_bending[gauss][i] * GPW #GPW = WGP * L     

        # coordinates
        x_1 = self.geometry.nodes[0].coordinates[0]
        y_1 = self.geometry.nodes[0].coordinates[1]
        x_2 = self.geometry.nodes[1].coordinates[0]
        y_2 = self.geometry.nodes[1].coordinates[1]

        # length, sin, cos
        length_0 = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2)   # original length
        c = (x_2 - x_1)/ length_0 # cos_alpha
        s = (y_2 - y_1)/ length_0 # sin_alpha

        # transformation matrix
        T = array([ [c, s, 0, 0, 0, 0],
                    [-s, c, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, c, s, 0],
                    [0, 0, 0, -s, c, 0],
                    [0, 0, 0, 0, 0, 1] ]) 

        RHS_glob = dot(transpose(T),RHS)
        print(RHS)
        return [LHS, RHS_glob]

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_X))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_Y))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[0], ROTATION))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[1], DISPLACEMENT_X))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[1], DISPLACEMENT_Y))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[1], ROTATION))  ## added by Andreas Riedl

        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_X))    ## added by Andreas Riedl
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_Y))    ## added by Andreas Riedl
        equation_ids.append(self.geometry[0].EquationId(ROTATION))    ## added by Andreas Riedl
        equation_ids.append(self.geometry[1].EquationId(DISPLACEMENT_X))    ## added by Andreas Riedl
        equation_ids.append(self.geometry[1].EquationId(DISPLACEMENT_Y))    ## added by Andreas Riedl
        equation_ids.append(self.geometry[1].EquationId(ROTATION))    ## added by Andreas Riedl

        return equation_ids

    def GetValues(self, step=0):
        values = zeros(self.GetDofsPerNode()*self.geometry.GetNumberOfNodes())
        values[0] = self.geometry[0].GetSolutionStepValue(DISPLACEMENT_X, step) ## added by Andreas Riedl
        values[1] = self.geometry[0].GetSolutionStepValue(DISPLACEMENT_Y, step) ## added by Andreas Riedl
        values[2] = self.geometry[0].GetSolutionStepValue(ROTATION, step) ## added by Andreas Riedl
        values[3] = self.geometry[1].GetSolutionStepValue(DISPLACEMENT_X, step) ## added by Andreas Riedl
        values[4] = self.geometry[1].GetSolutionStepValue(DISPLACEMENT_Y, step) ## added by Andreas Riedl
        values[5] = self.geometry[1].GetSolutionStepValue(ROTATION, step) ## added by Andreas Riedl

        return values
