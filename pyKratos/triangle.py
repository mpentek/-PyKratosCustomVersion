from __future__ import print_function, absolute_import, division 
import math
from numpy import *


class Triangle:

    def __init__(self, node_list):
        if(len(node_list) != 3):
            raise Exception("wrong number of nodes! should be 3!!")
        self.nodes = node_list

        for node in self.nodes:
            if(node.Id < 0):
                raise Exception("node with Id lesser than 0 found")

    # def Nodes(self):
        # return self.nodes

    def __getitem__(self, key):
        return self.nodes[key]
    
    def GetNumberOfNodes(self):
        return 3

    def GaussPoints(self, order):                   ## added by Andreas Riedl
        gpc = []
        weights = []

        x10 = self.nodes[1].coordinates[0] - self.nodes[0].coordinates[0]
        y10 = self.nodes[1].coordinates[1] - self.nodes[0].coordinates[1]
        x20 = self.nodes[2].coordinates[0] - self.nodes[0].coordinates[0]
        y20 = self.nodes[2].coordinates[1] - self.nodes[0].coordinates[1]

        detJ = x10 * y20 - y10 * x20

        if(order == 1):  # give back 1 single integration point
            gpc.append( [1/3, 1/3] ) #[xi, eta]
            Area = 0.5 * detJ
            weights = [Area]

        elif(order == 2):  # gives back 3 integration points
            gpc.append( [1/6, 1/6] )
            gpc.append( [2/3, 1/6] )
            gpc.append( [1/6, 2/3] )
            
            weights = [1/6 * detJ, 1/6 * detJ, 1/6 * detJ]
        
        elif(order == 3):  # gives back 4 integration points
            gpc.append( [1/3, 1/3] )
            gpc.append( [1/5, 1/5] )
            gpc.append( [1/5, 3/5] )
            gpc.append( [3/5, 1/5] )
            
            weights = [-27/48 * detJ, 25/48 * detJ, 25/48 * detJ, 25/48 * detJ]

        else:
            raise Exception("integration order not implemented")

        return [gpc, weights ]


    def ShapeFunctions(self, points):
        '''this function provides the shape function values, derivatives and integration_weight'''
        '''at the location of the gauss points. Order of integration is controlled'''
        '''by the optional parameter "order".'''
        '''N[gauss][i] contains the shape function of node i computed at the position of "gauss" '''
        '''derivatives[gauss][i,k] contains the derivative of node i, component k at the position of gauss '''
        '''weights[gauss] includes the integration weights, including the det of the jacobian, to be used '''
        '''at the gauss point'''
        derivatives = []
        Ncontainer = []

        x10 = self.nodes[1].coordinates[0] - self.nodes[0].coordinates[0]
        y10 = self.nodes[1].coordinates[1] - self.nodes[0].coordinates[1]
        x20 = self.nodes[2].coordinates[0] - self.nodes[0].coordinates[0]
        y20 = self.nodes[2].coordinates[1] - self.nodes[0].coordinates[1]

        detJ = x10 * y20 - y10 * x20 # =2 * Area

        DN_DX = zeros((3, 2), dtype=float)
        DN_DX[0, 0] = -y20 + y10    # -(y2-y0) + y1 - y0 -> b1
        DN_DX[0, 1] = x20 - x10     # x2 - x0 - (x1 -x0) -> c1
        DN_DX[1, 0] = y20           # y2 - y0   -> b2
        DN_DX[1, 1] = -x20          # -(x2 -x0) -> c2
        DN_DX[2, 0] = -y10          #-(y1-y0)   -> b3
        DN_DX[2, 1] = x10           #x1-x0  -> c3

        DN_DX /= detJ
        
        for point_i in points:
            Ncontainer.append( array([ 1 - point_i[0] - point_i[1], point_i[0], point_i[1] ]) )
            derivatives.append( DN_DX )

        return [Ncontainer, derivatives]
