from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = point2d.Point2D(list_of_nodes)
    return Point2DCondition(Id, prop, geom)

class Point2DCondition:

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry

    def GetDofsPerNode(self):
        return 3

                
    def CalculateLocalSystem(self,ProcessInfo, model_part):
        fx = self.geometry[0].GetSolutionStepValue(EXTERNAL_FORCE_X,0) 
        fy = self.geometry[0].GetSolutionStepValue(EXTERNAL_FORCE_Y,0)   
        moment = self.geometry[0].GetSolutionStepValue(EXTERNAL_MOMENT,0) 
        
        RHS = zeros(3)  
        LHS = zeros((3,3))
        
        RHS[0] = fx
        RHS[1] = fy
        RHS[2] = moment
        
        return [LHS, RHS]

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_X))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_Y))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[0], ROTATION))  ## added by Andreas Riedl
        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_X))    ## added by Andreas Riedl
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_Y))    ## added by Andreas Riedl
        equation_ids.append(self.geometry[0].EquationId(ROTATION))    ## added by Andreas Riedl
        return equation_ids

    def GetValues(self, step=0):
        values = zeros(self.GetDofsPerNode()*self.geometry.GetNumberOfNodes())
        values[0] = self.geometry[0].GetSolutionStepValue(DISPLACEMENT_X, step) ## added by Andreas Riedl
        values[1] = self.geometry[0].GetSolutionStepValue(DISPLACEMENT_Y, step) ## added by Andreas Riedl
        values[2] = self.geometry[0].GetSolutionStepValue(ROTATION, step) ## added by Andreas Riedl
        return values


