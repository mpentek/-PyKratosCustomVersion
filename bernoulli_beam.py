##                   pyKratos
##
##  License:		 BSD License
##					 pyKratos default license: pyKratos/LICENSE
##
##  Main authors:    Andreas Riedl, Technical University of Munich, Chair of Structural Analysis
##
##  contributors:    Máté Pentek, M.Sc., Technical University of Munich, Chair of Structural Analysis
##
##  description:     bernoulli beam element for the linear calculation of the deformation

from __future__ import print_function, absolute_import, division 
import math
from pyKratos import *
from numpy import *
from scipy import linalg

def Create(Id, prop, list_of_nodes):
    geom = line2d.Line2D(list_of_nodes)   
    return BernoulliBeamElement(Id, prop, geom)


class BernoulliBeamElement:
    integration_order = 3  # this is like a c "static variable" one for all of the objects of this type

    def __init__(self, Id, prop, geometry):
        self.Id = Id
        self.prop = prop
        self.geometry = geometry
        self.K_static_last = 0      # last calculated stiffness matrix
        self.T = 0                  # transformation matrix
        self.RHSstatic = 0          # load vector
        self.length = 0

        self.variables = []
        
        #mark as lagrangian all of the nodes in the element
        for node in self.geometry:
            node.SetSolutionStepValue(IS_LAGRANGIAN,0,True)

    def GetDofsPerNode(self):
        return 3
    
        #this _auxiliary function computes the stiffness contribution from bending
    def _ComputeStiffnessContributionBending(self,ProcessInfo):
        order = self.integration_order
        [gpc, weights] = self.geometry.GaussPoints(order)
        [Ns, dderivatives, ddderivatives] = self.geometry.ShapeFunctionsBending(gpc)        

        RHSstatic = zeros(4)    # preallocate load vector
        K = zeros((4, 4))       # preallocate stiffness matrix
               
        nnodes = self.geometry.GetNumberOfNodes()
        number_of_gauss = len(Ns)
        E = self.prop[YOUNG_MODULUS]    # in N/m^2
        I = self.prop[MOMENT_INERTIA_AREA]    # in m^4     

        # # preforming, imposed values - TO DO
        # u0_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_X"]
        # w0_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_Y"]
        # u0_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_X"]
        # w0_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_Y"]           

        #loop over all gauss points                
        for gauss in range(0, number_of_gauss):
            N = Ns[gauss]
            N_dderivate = dderivatives[gauss] #second derivate of N
            GPW = weights[gauss]

            # local stiffness matrix bending
            for i in range(0,4):
                for j in range(0,4):
                    K[i,j] += N_dderivate[i] * N_dderivate[j] * E*I * GPW   #GPW = WGP * L                  
        
        return [K, RHSstatic]

     #this _auxiliary function computes the stiffness contribution from stretching
    def _ComputeStiffnessContributionLinear(self,ProcessInfo):
        order = self.integration_order
        [gpc, weights] = self.geometry.GaussPoints(order=1)
        [Ns, derivatives] = self.geometry.ShapeFunctionsLinear(gpc)

        RHSstatic = zeros(2)    # preallocate load vector
        K = zeros((2, 2))       # preallocate stiffness matrix
               
        nnodes = self.geometry.GetNumberOfNodes()
        number_of_gauss = len(Ns)
        E = self.prop[YOUNG_MODULUS]    # in N/m^2
        A = self.prop[SECTION_TYPE]     # in m^        

        # # preforming, imposed values - TO DO
        # u0_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_X"]
        # w0_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_Y"]
        # u0_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_X"]
        # w0_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_Y"]

        #loop over all gauss points                
        for gauss in range(0, number_of_gauss):
            N = Ns[gauss]
            N_derivate = derivatives[gauss]
            GPW = weights[gauss]

            # local stiffness matrix
            for i in range(0,2):
                for j in range(0,2):
                    K[i,j] += N_derivate[i] * N_derivate[j] * E*A * GPW   #GPW = WGP * L               

        return [K, RHSstatic]

    
    # this function computes the selfweight of the beam and produce a line load condition
    def _CalculateSelfWeight(self, ProcessInfo, model_part):
        try:
            force_y = self.prop[BODY_FORCE_Y] #in m/s^2
            density = self.prop[DENSITY]      #in kg/m^3
            A = self.prop[SECTION_TYPE]     # in m^

            # coordinates
            x_1 = self.geometry.nodes[0].coordinates[0]
            y_1 = self.geometry.nodes[0].coordinates[1]
            x_2 = self.geometry.nodes[1].coordinates[0]
            y_2 = self.geometry.nodes[1].coordinates[1]

            length_0 = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2)   # original length
            F_Y_glob = density * A * length_0 * force_y

        except:
            F_Y_glob = 0
            print('SelfWeight failure')

        if F_Y_glob != 0:            
            
            #calculates lokal forces
            F_lok = dot(self.T[0:2,0:2],[0,F_Y_glob])
            F_X_lok = F_lok[0]
            F_Y_lok = F_lok[1]
            
            #set lokal line loads
            val_x = self.GetSolutionStepValue('LINE_LOAD_X',0)
            self.SetSolutionStepValue('LINE_LOAD_X', 0, F_X_lok + val_x)

            val_y = self.GetSolutionStepValue('LINE_LOAD_Y',0)
            self.SetSolutionStepValue('LINE_LOAD_Y', 0, F_Y_lok + val_y)

            # creates new line_load_condition
            if val_x == 0 and val_y == 0: #if not already exist
                try: new_cond_Id = max(model_part.Conditions.keys())+1
                except: new_cond_Id = 0
                cond = { new_cond_Id: [0,[self.Id]] }
                model_part.AddElementConditions('line_load_condition_2d',cond)
            

    # called from each element and added to the system matrix
    def CalculateLocalSystem(self,ProcessInfo, model_part):
        self._CalculateTransformationMatrix()

        self._CalculateSelfWeight(ProcessInfo, model_part)

        K_beam, RHSstatic_beam = self._ComputeStiffnessContributionBending(ProcessInfo)
        K_truss, RHSstatic_truss = self._ComputeStiffnessContributionLinear(ProcessInfo)

        # build 6x6 local stiffness matrix
        K = K_beam
        K = insert(K,[0,2],0,1) #from 4x4 to 6x6
        K = insert(K,[0,2],0,0)
        K[0,0] = K_truss[0,0]
        K[3,0] = K_truss[1,0]
        K[0,3] = K_truss[0,1]
        K[3,3] = K_truss[1,1]

        # build 6x1 local rhs
        RHSstatic = RHSstatic_beam
        RHSstatic = insert(RHSstatic,0,RHSstatic_truss[0],0)
        RHSstatic = insert(RHSstatic,3,RHSstatic_truss[1],0)

        # global stiffness matrix
        K_glob = dot(transpose(self.T),dot(K,self.T))
        
        # global load vector
        RHSstatic_glob = dot(transpose(self.T),RHSstatic)

        # preforming
        # preforming = array([u0_1, w0_1, u0_2, w0_2])   
        # RHS_preforming = -dot(K_glob, preforming)        
        # RHSstatic_glob = RHSstatic_glob + RHS_preforming

        self.RHSstatic = RHSstatic # save newest local rhs for sgr calculation
        self.K_static_last = K # save newest local K for sgr calculation

        return [K_glob,RHSstatic_glob]


        # this function calculates the transformation matrix from local to global
    def _CalculateTransformationMatrix(self):
        # coordinates
        x_1 = self.geometry.nodes[0].coordinates[0]
        y_1 = self.geometry.nodes[0].coordinates[1]
        x_2 = self.geometry.nodes[1].coordinates[0]
        y_2 = self.geometry.nodes[1].coordinates[1]

        # length, sin, cos
        length_0 = math.sqrt((x_2 - x_1)**2 + (y_2 - y_1)**2)   # original length
        self.length = length_0   
        c = (x_2 - x_1)/ length_0 # cos_alpha
        s = (y_2 - y_1)/ length_0 # sin_alpha
        alpha = (arcsin(s)*180)/pi

        # transformation matrix
        T = array([ [c, s, 0, 0, 0, 0],
                    [-s, c, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, c, s, 0],
                    [0, 0, 0, -s, c, 0],
                    [0, 0, 0, 0, 0, 1] ]) 

        self.T = T


        # this function returns a list with the node and unkowns to be solved for
    def GetDofList(self):
        unknowns = []
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_X))  
        unknowns.append(Dof(self.geometry[0], DISPLACEMENT_Y)) 
        unknowns.append(Dof(self.geometry[0], ROTATION)) 
        unknowns.append(Dof(self.geometry[1], DISPLACEMENT_X))  
        unknowns.append(Dof(self.geometry[1], DISPLACEMENT_Y))  
        unknowns.append(Dof(self.geometry[1], ROTATION)) 
        return unknowns

        # this function returns a list with the equationIds to be solved for
    def EquationId(self):   
        equation_ids = []
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_X))
        equation_ids.append(self.geometry[0].EquationId(DISPLACEMENT_Y))  
        equation_ids.append(self.geometry[0].EquationId(ROTATION))  
        equation_ids.append(self.geometry[1].EquationId(DISPLACEMENT_X))    
        equation_ids.append(self.geometry[1].EquationId(DISPLACEMENT_Y))  
        equation_ids.append(self.geometry[1].EquationId(ROTATION))    
        return equation_ids


    #calculates the local displacments from the global
    def CalculateLocalDisplacement(self):
        # displacement
        u_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_X"]
        w_1 = self.geometry.nodes[0].variables[0]["DISPLACEMENT_Y"]
        r_1 = self.geometry.nodes[0].variables[0]["ROTATION"]
        u_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_X"]
        w_2 = self.geometry.nodes[1].variables[0]["DISPLACEMENT_Y"]
        r_2 = self.geometry.nodes[1].variables[0]["ROTATION"]

        # transform to lokal
        u_glob = [u_1, w_1, r_1, u_2, w_2, r_2]
        u_lok = dot(self.T,u_glob)

        # internal forces at the ends
        s_lok = dot(self.K_static_last,u_lok) - self.RHSstatic
        sgr_lok = s_lok
        sgr_lok[0]=sgr_lok[0]*-1
        sgr_lok[1]=sgr_lok[1]*-1
        sgr_lok[2]=sgr_lok[2]*-1

        return [u_lok]

    # Calculates the internal forces N,V,M for 50 steps of the beam
    def CalculateInternalForces(self):

        [u_lok] = self.CalculateLocalDisplacement() #calculate local displacement
        xi = linspace(0,1,num=51)
        x = xi * self.length

        N = zeros(51)
        V = zeros(51)
        M = zeros(51)

        E = self.prop[YOUNG_MODULUS]    # in N/m^2
        A = self.prop[SECTION_TYPE]     # in m^    
        I = self.prop[MOMENT_INERTIA_AREA]    # in m^4  

        [Ns_linear, derivatives_linear] = self.geometry.ShapeFunctionsLinear(xi)
        [Ns_bend, dderivatives_bend, ddderivatives_bend] = self.geometry.ShapeFunctionsBending(xi)

        #internal forces from Boundaries
        for i in range(0,51):
            M[i] = -1 * E * I *( u_lok[1]*dderivatives_bend[i][0] + u_lok[2]*dderivatives_bend[i][1] + u_lok[4]*dderivatives_bend[i][2] + u_lok[5]*dderivatives_bend[i][3] )
            V[i] = -1 * E * I *( u_lok[1]*ddderivatives_bend[i][0] + u_lok[2]*ddderivatives_bend[i][1] + u_lok[4]*ddderivatives_bend[i][2] + u_lok[5]*ddderivatives_bend[i][3] )
            N[i] = E * A * (u_lok[0]*derivatives_linear[i][0] + u_lok[3]*derivatives_linear[i][1])

            qx = self.GetSolutionStepValue(LINE_LOAD_X,0)
            qy = self.GetSolutionStepValue(LINE_LOAD_Y,0)  

        #forces from const loads 
        for i in range(0,51):
            M[i] += -1 * E * I * qy * self.length**2 / (24 * E * I) * (2 - 12 *xi[i] + 12 * xi[i]**2)
            V[i] += ((-1) * E * I * qy * self.length**2 / (24 * E * I)) * (1 / self.length) * (-12 + 24 * xi[i])
            N[i] += E * A * qx * self.length / (E * A * 2) * (1 - 2 * xi[i])       
   
        return [N,V,M]


    # calculates the coordinates of the deformed beam
    def CalculateDeflectionCurve(self, scale):

        [u_lok] = self.CalculateLocalDisplacement() #calculate local displacement
        xi = linspace(0, 1, num=51)

        E = self.prop[YOUNG_MODULUS]    # in N/m^2
        A = self.prop[SECTION_TYPE]     # in m^    
        I = self.prop[MOMENT_INERTIA_AREA]    # in m^4  

        u_curve_local = zeros(51)
        w_curve_local = zeros(51)

        [Ns_linear, derivatives_linear] = self.geometry.ShapeFunctionsLinear(xi)
        [Ns_bend, dderivatives_bend, ddderivatives_bend] = self.geometry.ShapeFunctionsBending(xi)

        qx = self.GetSolutionStepValue(LINE_LOAD_X,0)
        qy = self.GetSolutionStepValue(LINE_LOAD_Y,0)  

        for i in range(0,51):
            # deformation from boundaries
            u_curve_local[i] = u_lok[0]*Ns_linear[i][0] + u_lok[3]*Ns_linear[i][1]
            w_curve_local[i] = u_lok[1]*Ns_bend[i][0] + u_lok[2]*Ns_bend[i][1] + u_lok[4]*Ns_bend[i][2] + u_lok[5]*Ns_bend[i][3]

            # deformation from const load
            u_curve_local[i] += qx * self.length**2 / (E * A * 2) * ( xi[i]- xi[i]**2)
            w_curve_local[i] += qy * self.length**2 * self.length**2 / (24 * E * I) * xi[i]**2 * (1 - 2 * xi[i] + xi[i]**2)

        # transform to global
        curve_global = dot(transpose(self.T[0:2,0:2]), array((u_curve_local,w_curve_local)))
        u_curve_global = curve_global[0]
        w_curve_global = curve_global[1]

        # coordinates for global plot - x_def and y_def
        x_1 = self.geometry.nodes[0].coordinates[0]
        y_1 = self.geometry.nodes[0].coordinates[1]
        x_2 = self.geometry.nodes[1].coordinates[0]
        y_2 = self.geometry.nodes[1].coordinates[1]

        x_def = zeros(51)
        y_def = zeros(51)
        for i in range(0,51):
            x_def[i] = x_1 * (1-xi[i]) + x_2 * xi[i] + u_curve_global[i] * scale
            y_def[i] = y_1 * (1-xi[i]) + y_2 * xi[i] + w_curve_global[i] * scale

        return [x_def, y_def]
    

    def GetSolutionStepValue(self, variable_name, step):
        try: return self.variables[step][variable_name]
        except: return 0 

    def SetSolutionStepValue(self, variable_name, step, value):        
        try: self.variables[step][variable_name] = value
        except:
            self.variables.append(dict())
            self.variables[step][variable_name] = value
