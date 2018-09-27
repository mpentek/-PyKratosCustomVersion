from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create( Id, prop, list_of_nodes):
    geom = triangle.Triangle(list_of_nodes)
    return SurfaceLoad2DCondition(Id, prop, geom)

class SurfaceLoad2DCondition:

    def __init__(self, Id, prop, geom):
        self.Id = Id
        self.prop = prop
        self.geometry = geom
        self.element = []

    def GetDofsPerNode(self):
        return 2

                
    def CalculateLocalSystem(self, ProcessInfo, model_part):
        list_of_nodes = self.geometry.nodes

        for ele in model_part.Elements.values():
            if ele.geometry.nodes == list_of_nodes:
                self.element = ele   

        try: px = self.element.GetSolutionStepValue(SURFACE_LOAD_X,0)
        except: px = 0
        try: py = self.element.GetSolutionStepValue(SURFACE_LOAD_Y,0)  
        except: py=0 

        # for constant surface loads. otherwise surface load input must be for each node
        p_x_load = array([ px, px, px ])     #p_x_1, p_x_2, p_x_3
        p_y_load = array([ py, py, py ])     #p_y_1, p_y_2, p_y_3
                
        RHS = zeros(6)  
        LHS = zeros((6,6))

        order = 2
        [gpc, weights] = self.geometry.GaussPoints(order)
        [Ns, derivatives] = self.geometry.ShapeFunctions(gpc)

        number_of_gauss = len(gpc)     
        
         #loop over all gauss points                
        for gauss in range(0, number_of_gauss):
            N = Ns[gauss]

            # computation of external surface load
            p_x = dot(N, p_x_load)
            p_y = dot(N, p_y_load)

            RHS_x = weights[gauss] * dot( N.transpose(), p_x)
            RHS_y = weights[gauss] * dot( N.transpose(), p_y)
            RHS[0:3] += RHS_x
            RHS[3:6] += RHS_y   

        return [LHS, RHS]

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_X))
        unknowns.append(Dof(self.geometry[1], DISPLACEMENT_X))
        unknowns.append(Dof(self.geometry[2], DISPLACEMENT_X))
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_Y))
        unknowns.append(Dof(self.geometry[1], DISPLACEMENT_Y))
        unknowns.append(Dof(self.geometry[2], DISPLACEMENT_Y))

        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_X))
        equation_ids.append(self.geometry[1].EquationId(DISPLACEMENT_X))
        equation_ids.append(self.geometry[2].EquationId(DISPLACEMENT_X))
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_Y))
        equation_ids.append(self.geometry[1].EquationId(DISPLACEMENT_Y))
        equation_ids.append(self.geometry[2].EquationId(DISPLACEMENT_Y))

        return equation_ids

    def GetValues(self, step=0):
        values = zeros(self.GetDofsPerNode()*self.geometry.GetNumberOfNodes())
        values[0] = self.geometry[0].GetSolutionStepValue(DISPLACEMENT_X, step) ## added by Andreas Riedl
        values[1] = self.geometry[1].GetSolutionStepValue(DISPLACEMENT_X, step) ## added by Andreas Riedl
        values[2] = self.geometry[2].GetSolutionStepValue(DISPLACEMENT_X, step) ## added by Andreas Riedl
        values[3] = self.geometry[3].GetSolutionStepValue(DISPLACEMENT_Y, step) ## added by Andreas Riedl
        values[4] = self.geometry[4].GetSolutionStepValue(DISPLACEMENT_Y, step) ## added by Andreas Riedl
        values[5] = self.geometry[5].GetSolutionStepValue(DISPLACEMENT_Y, step) ## added by Andreas Riedl

        return values


