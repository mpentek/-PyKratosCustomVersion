##                   pyKratos
##
##  License:		 BSD License
##					 pyKratos default license: pyKratos/LICENSE
##
##  Main authors:    Andreas Riedl, Technical University of Munich, Chair of Structural Analysis
##
##  contributors:    Máté Pentek, M.Sc., Technical University of Munich, Chair of Structural Analysis
##
##  description:     truss element for the geometrical nonlinear calculation with selfweight

from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg


def Create(Id, prop, list_of_nodes):
    geom = line2d.Line2D(list_of_nodes)   
    return TrussElement(Id, prop, geom)


class TrussElement:
    integration_order = 2  # this is like a c "static variable" one for all of the objects of this type

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry
        self.K_static_last = 0
        self.T = 0
       
        #mark as lagrangian all of the nodes in the element
        for node in self.geometry:
            node.SetSolutionStepValue(IS_LAGRANGIAN,0,True)

    def GetDofsPerNode(self):
        return 2
    
        #this _auxiliary function computes the stiffness contribution
    def _ComputeStiffnessContribution(self,ProcessInfo):
        order = self.integration_order
        [gpc, weights] = self.geometry.GaussPoints(order)
        [Ns, derivatives] = self.geometry.ShapeFunctionsLinear(gpc)

        RHSstatic = zeros(4)    # preallocate load vector
        K_ele = zeros((4, 4))   # preallocate stiffness matrix elastic part
        K_geo = zeros((4, 4))   # preallocate stiffness matrix geometric part
               
        nnodes = self.geometry.GetNumberOfNodes()
        number_of_gauss = len(Ns)
        E = self.prop[YOUNG_MODULUS]    # in N/m^2
        A = self.prop[SECTION_TYPE]     # in m^2

        # coordinates
        x_1 = self.geometry.nodes[0].coordinates[0]
        y_1 = self.geometry.nodes[0].coordinates[1]
        x_2 = self.geometry.nodes[1].coordinates[0]
        y_2 = self.geometry.nodes[1].coordinates[1]

        # displacement
        u_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_X"]
        w_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_Y"]
        u_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_X"]
        w_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_Y"]

        # preforming, imposed values    NOT IMPLEMENTED!
        #u0_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_X"]
        #w0_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_Y"]
        #u0_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_X"]
        #w0_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_Y"]

        # length, sin, cos
        length_0 = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2)   # original length
        length = math.sqrt((x_2+u_2-x_1-u_1)**2 + (y_2+w_2-y_1-w_1)**2) # current length
        #print([length_0, length])
        c = (x_2 + u_2 - x_1 - u_1)/ length # cos_alpha
        s = (y_2 + w_2 - y_1 - w_1)/ length # sin_alpha    
        alpha = (arcsin(s)*180)/pi 


        # transformation matrix
        T = array([ [c, s, 0, 0],
                    [-s, c, 0, 0],
                    [0, 0, c, s],
                    [0, 0, -s, c] ])

        self.T = T

        tension = (E*A)/length_0 *((u_2-u_1)*c + (w_2-w_1)*s) + (E*A)/(2*length_0**2) * ((u_1-u_2)*s + (w_2-w_1)*c)**2

        #loop over all gauss points
        for gauss in range(0, number_of_gauss):
            N = Ns[gauss]
            N_derivate = derivatives[gauss]
            GPW = weights[gauss]

            # local stiffness matrix - elastic part
            for i in range(0,2):
                for j in range(0,2):
                    K_ele[2*i,2*j] += N_derivate[i] * N_derivate[j] * E*A * GPW   #GPW = WGP * L            
            
            # local stiffness matrix - geometric part
            for i in range(0,2):
                for j in range(0,2):
                    K_geo[(2*i)+1,(2*j)+1] += N_derivate[i] * N_derivate[j] * tension * GPW * length_0/length   #GPW = WGP * L           
     
            # local load vector - qx TRANSFER NOT IMPLEMENTED
            qx = 0
            for i in range(0, 2):
                RHSstatic[2*i] += qx * N[i] * GPW #GPW = WGP * L            
        
        # local stiffness matrix
        K_static = K_ele + K_geo
        self.K_static_last = K_static # save newest K for sgr calculation
        # global stiffness matrix
        K_glob = dot(transpose(T),dot(K_static,T))
        
        # global load vector
        RHSstatic_glob = dot(transpose(T),RHSstatic)

        ## preforming - NOT IMPLEMENTED
        #preforming = array([u0_1, w0_1, u0_2, w0_2])        
        #RHS_preforming = dot(K_glob, preforming)        
        #RHSstatic_glob = RHSstatic_glob + RHS_preforming
        
        return [K_glob, RHSstatic_glob]  
    
        # called from each element and added to the system matrix
    def CalculateLocalSystem(self,ProcessInfo, model_part):
        K_glob, RHSstatic = self._ComputeStiffnessContribution(ProcessInfo)
               
        RHS = RHSstatic
        LHS = K_glob

        
        return [LHS,RHS]

        # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_X))  
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_Y))  
        unknowns.append(Dof(self.geometry[1], DISPLACEMENT_X))  
        unknowns.append(Dof(self.geometry[1], DISPLACEMENT_Y))  
        return unknowns

        # this function returns a list with the equationIds to be solved for
    def EquationId(self):   
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_X))
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_Y))   
        equation_ids.append(self.geometry[1].EquationId(DISPLACEMENT_X))    
        equation_ids.append(self.geometry[1].EquationId(DISPLACEMENT_Y))    
        return equation_ids

    # calculates the internal forces of the truss - only N
    def CalculateInternalForces(self):

        # displacement
        u_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_X"]
        w_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_Y"]
        u_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_X"]
        w_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_Y"]

        # internal forces
        u_glob = [u_1, w_1, u_2, w_2]
        u_lok = dot(self.T,u_glob)
        s_lok = dot(self.K_static_last, u_lok)

        sgr_lok = s_lok
        sgr_lok[0]=sgr_lok[0]*-1
        sgr_lok[1]=sgr_lok[1]*-1
        N = (sgr_lok[0]+sgr_lok[2])/2   

        return N

    def CalculateDeflectionCurve(self, scale):

        x_1 = self.geometry.nodes[0].coordinates[0]
        y_1 = self.geometry.nodes[0].coordinates[1]
        x_2 = self.geometry.nodes[1].coordinates[0]
        y_2 = self.geometry.nodes[1].coordinates[1]

        # displacement
        u_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_X"]
        w_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_Y"]
        u_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_X"]
        w_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_Y"]

        # new global coordinates
        x_1_new = x_1 + u_1 * scale
        y_1_new = y_1 + w_1 * scale
        x_2_new = x_2 + u_2 * scale
        y_2_new = y_2 + w_2 * scale

        x_def = array((x_1_new,x_2_new))
        y_def = array((y_1_new,y_2_new))        

        return [x_def, y_def]

