from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = triangle.Triangle(list_of_nodes)
    return TriangularShell2D(Id, prop, geom)


class TriangularShell2D:
    integration_order = 2  # this is like a c "static variable" one for all of the objects of this type

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry

        self.variables = []

        
        #mark as lagrangian all of the nodes in the element
        for node in self.geometry:
            node.SetSolutionStepValue(IS_LAGRANGIAN,0,True)

    def GetDofsPerNode(self):
        return 2



    

    #this _auxiliary function computes the stiffness contribution
    def _ComputeStiffnessContribution(self,ProcessInfo):
        order = self.integration_order
        [gpc, weights] = self.geometry.GaussPoints(order)
        [Ns, derivatives] = self.geometry.ShapeFunctions(gpc)

        RHSstatic = zeros(6)  # no external forces so far
        K = zeros((6, 6))
                
        nnodes = self.geometry.GetNumberOfNodes()

        number_of_gauss = len(Ns)
        for gauss in range(0, number_of_gauss):
            N = Ns[gauss]
            DN_DX = derivatives[gauss]

            # computation of internal forces
            C = self.ComputeElasticityTensor() # material tensor
            B = self.ComputeB(DN_DX)    # ableitungstensor

            height = 0.2    # height of shell element   ##TO DO
            tmp = dot(C, B)
            K += height * weights[gauss] * dot(B.transpose(), tmp) 

           

            # # computation of external line loads
            # q=[]
            # q_1 = array([ 1, 1, 1, 1 ])   #q12_x, q13_x, q12_y, q13_y
            # q_2 = array([ 1, 1, 1, 1 ])   #q21_x, q23_x, q21_y, q23_y
            # q_3 = array([ 1, 1, 1, 1 ])   #q31_x, q32_x, q31_y, q32_y

            # N1 = [ [N[2], N[3], 0, 0], [0, 0, N[2], N[3]] ]
            # N2 = [ [N[1], N[3], 0, 0], [0, 0, N[1], N[3]] ]
            # N3 = [ [N[1], N[2], 0, 0], [0, 0, N[1], N[2]] ]

            # q_1 = dot( N1, q_1 )
            # q_2 = dot( N2, q_2 )
            # q_3 = dot( N3, q_3 )

            # N1 = [0, [N[2], N[3], 0, 0, 0], [0, 0, 0, 0, N[2], N[3]] ]
            # N2 = [ [N[1], 0, N[3], 0, 0, 0], [0, 0, 0, N[1], 0, N[3]] ]
            # N3 = [ [N[1], N[2], 0, 0, 0, 0], [0, 0, 0, N[1], N[2], 0] ]            

           
            

        # compute RESIDUAL for the RHS
        # since the problem is linear this can be done as
        # RHS = fext - LHS*values
        # note that this is done out of the integration loop!
        #values = self._GetAllValuesVector(DISPLACEMENT_X,DISPLACEMENT_Y)  # get values of unknown at the nodes
        #RHSstatic -= dot(K, values)

        return [K, RHSstatic]
        
    def CalculateLocalSystem(self,ProcessInfo, model_part):

        K, RHSstatic = self._ComputeStiffnessContribution(ProcessInfo)
        
        #get BDF coefficients to correctly compute the LHS
        # coeffs =  ProcessInfo[BDF_COEFFICIENTS]
        # c0 = coeffs[0]
        # c1 = coeffs[1]
        # c2 = coeffs[2]

        
        #RHS = RHSstatic + RHSinertia
        #LHS = c0*M + 1.0/c0*K
        LHS = K
        RHS = RHSstatic
        
        return [LHS,RHS]

    def ComputeB(self, DN_DX):

        B = zeros((3, 6))

        for i in range(0, 3):
            B[0, i] = DN_DX[i, 0]
            B[2, i] = DN_DX[i, 1]
            B[1, i+3] = DN_DX[i, 1]
            B[2, i+3] = DN_DX[i, 0]
           
        return B

    def ComputeElasticityTensor(self):
        C = zeros((3, 3))
        E = self.prop[YOUNG_MODULUS]
        NU = self.prop[POISSON_RATIO]

        c1 = E / (1.00 - NU * NU)

        c2 = NU * c1

        c3 = E / 2 / (1 + NU)

        C = [[c1, c2, 0.0], [c2, c1, 0.0], [0.0, 0.0, c3]]

        return C

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

    def _GetAllValuesVector(self, var_x,var_y, step=0):
        values = zeros(6)
        values[0] = self.geometry[0].GetSolutionStepValue(var_x, step)
        values[2] = self.geometry[1].GetSolutionStepValue(var_x, step)
        values[4] = self.geometry[2].GetSolutionStepValue(var_x, step)
        values[1] = self.geometry[0].GetSolutionStepValue(var_y, step)
        values[3] = self.geometry[1].GetSolutionStepValue(var_y, step)
        values[5] = self.geometry[2].GetSolutionStepValue(var_y, step)
        return values


    # def _GetAccelerationsVector(self,ProcessInfo):
    #     values = zeros(6)
    #     coeffs =  ProcessInfo[BDF_COEFFICIENTS]
    #     c0 = coeffs[0]
    #     c1 = coeffs[1]
    #     c2 = coeffs[2]
                
    #     v0 = self._GetAllValuesVector(VELOCITY_X,VELOCITY_Y,0)
    #     v1 = self._GetAllValuesVector(VELOCITY_X,VELOCITY_Y,1)
    #     v2 = self._GetAllValuesVector(VELOCITY_X,VELOCITY_Y,2)
               
    #     return c0*v0+c1*v1+c2*v2

    def GetSolutionStepValue(self, variable_name, step):
        try: return self.variables[step][variable_name]
        except: return 0 

    def SetSolutionStepValue(self, variable_name, step, value):        
        try: self.variables[step][variable_name] = value
        except:
            self.variables.append(dict())
            self.variables[step][variable_name] = value
