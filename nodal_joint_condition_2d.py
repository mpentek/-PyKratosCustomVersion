from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = line2d.Line2D(list_of_nodes)
    return Joint2DCondition(Id, prop, geom)

class Joint2DCondition:

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry

    def GetDofsPerNode(self):
        return 3

                
    def CalculateLocalSystem(self,ProcessInfo, model_part):
        LHS = []
        RHS = []
        
        return [LHS, RHS]

    # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_X))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_Y))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[0], ROTATION))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_X))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_Y))  ## added by Andreas Riedl
        unknowns.append(Dof(self.geometry[1], ROTATION))  ## added by Andreas Riedl
        return unknowns

    def EquationId(self):
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_X))    ## added by Andreas Riedl
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_Y))    ## added by Andreas Riedl
        return equation_ids

    def GetValues(self, step=0):
        values = zeros(self.GetDofsPerNode()*self.geometry.GetNumberOfNodes())
        values[0] = self.geometry[0].GetSolutionStepValue(DISPLACEMENT_X, step) ## added by Andreas Riedl
        values[1] = self.geometry[0].GetSolutionStepValue(DISPLACEMENT_Y, step) ## added by Andreas Riedl
        return values


