from __future__ import print_function, absolute_import, division 
import math
from numpy import *


class Line2D:

    def __init__(self, node_list):
        if(len(node_list) != 2):
            raise Exception("wrong number of nodes! should be 2!!")
        self.nodes = node_list

        for node in self.nodes:
            if(node.Id < 0):
                raise Exception("node with Id lesser than 0 found")

    # def Nodes(self):
        # return self.nodes

    def __getitem__(self, key):
        return self.nodes[key]
    
    def GetNumberOfNodes(self):
        return 2

    def ShapeFunctionsLinear(self, points):
        '''this function provides the shape function values, derivatives and integration_weight'''
        '''at the location of the gauss points. Order of integration is controlled'''
        '''by the optional parameter "order".'''
        '''N[gauss][i] contains the shape function of node i computed at the position of "gauss" '''
        '''derivatives[gauss][i,k] contains the derivative of node i, component k at the position of gauss '''
        '''weights[gauss] includes the integration weights, including the det of the jacobian, to be used '''
        '''at the gauss point'''
        derivatives = [] #derivative
        Ncontainer = []        

        x10 = self.nodes[1].coordinates[0] - self.nodes[0].coordinates[0]
        y10 = self.nodes[1].coordinates[1] - self.nodes[0].coordinates[1]

        length = math.sqrt(x10**2 + y10**2)

        for xi in points:                            
                Ncontainer.append(array([ 1-xi, xi]))
                derivatives.append(array([ -1/length, 1/length ]))
                
        return [Ncontainer, derivatives]

    def ShapeFunctionsBending(self, points):                         ## added by Andreas Riedl
        '''this function provides the bending shape function values, second derivatives, third derivatives and integration_weight'''
        '''at the location of the gauss points. Order of integration is controlled'''
        '''by the optional parameter "order".'''
        '''N[gauss][i] contains the shape function of node i computed at the position of "gauss" '''
        '''derivatives[gauss][i,k] contains the derivative of node i, component k at the position of gauss '''
        '''weights[gauss] includes the integration weights, including the det of the jacobian, to be used '''
        '''at the gauss point'''
        dderivatives = [] #second derivative
        ddderivatives = []  #third derivative
        Ncontainer = []      

        x10 = self.nodes[1].coordinates[0] - self.nodes[0].coordinates[0]
        y10 = self.nodes[1].coordinates[1] - self.nodes[0].coordinates[1]
        length = math.sqrt(x10**2 + y10**2)  

        for xi in points:                            
                Ncontainer.append(array([ 1 - 3 * xi**2 + 2 * xi**3, (-xi + 2 * xi**2 - xi**3) * length, 3 * xi**2 - 2 * xi**3, (xi**2 - xi**3) * length ]))
                dderivatives.append(array([ (-6 / length**2 + 12 * xi / length**2), (4 / length - 6 * xi / length), (6 / length**2 - 12 * xi / length**2), (2 / length - 6 * xi / length) ]))
                ddderivatives.append(array([ 12 / length**3, -6 / length**2, -12 / length**3, -6 / length**2 ]))

        return [Ncontainer, dderivatives, ddderivatives]

    def GaussPoints(self,order):                         ## added by Andreas Riedl
        gpc = []
        weights = []

        x10 = self.nodes[1].coordinates[0] - self.nodes[0].coordinates[0]
        y10 = self.nodes[1].coordinates[1] - self.nodes[0].coordinates[1]

        length = math.sqrt(x10**2 + y10**2)
        
        if(order == 1):  # give back 1 single integration point
            gpc.append(-math.sqrt(3.0 / 5.0) / 2.0 + 0.5)
            weights = [length]            

        elif(order == 2):  # gives back 2 integration points
            gpc.append(-math.sqrt(1.0 / 3.0) / 2.0 + 0.5)
            gpc.append(math.sqrt(1.0 / 3.0) / 2.0 + 0.5)

            weights = [0.5*length, 0.5*length]           

        elif(order == 3):  # gives back 3 integration points       
            gpc.append(-math.sqrt(3.0 / 5.0) / 2.0 + 0.5)
            gpc.append(0.5)
            gpc.append(math.sqrt(3.0 / 5.0) / 2.0 + 0.5)         

            weights = [(5/18)*length, (8/18)*length, (5/18)*length]   

        elif(order == 4):  # gives back 4 integration points
            gpc.append( 1 / 2 * (-(3 / 7 + 2 / 7 * (6 / 5) ** 0.5) ** 0.5 + 1))
            gpc.append( 1 / 2 * (-(3 / 7 - 2 / 7 * (6 / 5) ** 0.5) ** 0.5 + 1))
            gpc.append( 1 / 2 * ((3 / 7 - 2 / 7 * (6 / 5) ** 0.5) ** 0.5 + 1)) 
            gpc.append( 1 / 2 * ((3 / 7 + 2 / 7 * (6 / 5) ** 0.5) ** 0.5 + 1))        

            weights = [(1 / 2 * (18 - (30) ** 0.5) / 36)*length, (1 / 2 * (18 + (30) ** 0.5) / 36)*length, (1 / 2 * (18 + (30) ** 0.5) / 36)*length, (1 / 2 * (18 - (30) ** 0.5) / 36)*length]  

        elif(order == 5):  # gives back 5 integration points      
            gpc.append( 1 / 2 * (-1 / 3 * (5 + 2 * (10 / 7) ** 0.5) ** 0.5 + 1))
            gpc.append( 1 / 2 * (-1 / 3 * (5 - 2 * (10 / 7) ** 0.5) ** 0.5 + 1))
            gpc.append( 1 / 2 * (0 + 1)) 
            gpc.append( 1 / 2 * (1 / 3 * (5 - 2 * (10 / 7) ** 0.5) ** 0.5 + 1))        
            gpc.append( 1 / 2 * (1 / 3 * (5 + 2 * (10 / 7) ** 0.5) ** 0.5 + 1))

            weights = [(1 / 2 * (322 - 13 * (70) ** 0.5) / 900)*length, (1 / 2 * (322 + 13 * (70) ** 0.5) / 900)*length, (1 / 2 * 128 / 225)*length, (1 / 2 * (322 + 13 * (70) ** 0.5) / 900)*length, (1 / 2 * (322 - 13 * (70) ** 0.5) / 900)*length]  

        else:
            raise Exception("integration order not implemented")

        return [gpc, weights]





