##                   pyKratos
##
##  License:		 BSD License
##					 pyKratos default license: pyKratos/LICENSE
##
##  Main authors:    Andreas Riedl, Technical University of Munich, Chair of Structural Analysis
##
##  contributors:    Máté Pentek, M.Sc., Technical University of Munich, Chair of Structural Analysis


from __future__ import print_function, absolute_import, division 
import numpy as np
from scipy.interpolate import griddata
import numpy.ma as ma
from pyKratos import *

import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def PlotInternalForces(Elements):
    xi = linspace(0,1,num=51)

    fig = plt.figure()  
    nelements = len(Elements) 
    i=1
    sub = {}
    ax = {}

    n_scale = 0.001 #from N to kN
    v_scale = 0.001
    m_scale = 0.001

    for element in Elements:
        [N,V,M] = element.CalculateInternalForces()

        #normal force
        sub[i] = fig.add_subplot(nelements, 3, i)
        if i == 1:
            sub[i].set_title('normal force')  

        val = N*n_scale
        sub[i].plot(xi, val, 'k')
        sub[i].set_ylabel('element '+ str(element.Id))
        sub[i].grid()
        sub[i].set_xlim(-0.1, 1.1)
        y_max = max( max(val), -min(val) )
        sub[i].set_ylim(-y_max*1.1, y_max*1.1)

        plt.fill_between(xi, val, 0, color = 'gray')
        i+=1

        #shear force
        sub[i] = fig.add_subplot(nelements, 3, i)
        if i == 2:
            sub[i].set_title('shear force')

        val = V*v_scale
        sub[i].plot(xi, val, 'k')
        sub[i].grid()
        sub[i].set_xlim(-0.1, 1.1)
        y_max = max( max(val), -min(val) )
        sub[i].set_ylim(y_max*1.1, -y_max*1.1)
        plt.fill_between(xi, val, 0, color = 'gray')
        i+=1

        #bending moment
        sub[i] = fig.add_subplot(nelements, 3, i)
        if i == 3:
            sub[i].set_title('bending moment')

        val = M*m_scale
        sub[i].plot(xi, val, 'k')
        sub[i].grid()
        sub[i].set_xlim(-0.1, 1.1)
        y_max = max( max(val), -min(val) )
        sub[i].set_ylim(y_max*1.1, -y_max*1.1)
        plt.fill_between(xi, val, 0, color = 'gray')
        i+=1

    plt.show()


    ###### systemplots with forces

    # creates 3 plots
    fig2 = plt.figure()
    sub[0] = fig2.add_subplot(1, 3, 1)
    sub[0].set_title('normal force')
    sub[1] = fig2.add_subplot(1, 3, 2)
    sub[1].set_title('shear force')
    sub[2] = fig2.add_subplot(1, 3, 3)
    sub[2].set_title('moment force')
    val_local ={}

    x = []  # preallocation
    y = []

    #over each element
    for element in Elements:
        [N,V,M] = element.CalculateInternalForces()
        val_local[0] = N * n_scale *0.1
        val_local[1] = V * v_scale *0.1
        val_local[2] = M * m_scale *0.1

        # coordinates for global plot - x_def and y_def
        x_1 = element.geometry.nodes[0].coordinates[0]
        y_1 = element.geometry.nodes[0].coordinates[1]
        x_2 = element.geometry.nodes[1].coordinates[0]
        y_2 = element.geometry.nodes[1].coordinates[1]

        x_org = linspace(x_1, x_2, num=51) 
        y_org = linspace(y_1, y_2, num=51)      

        # calc and plot n, v and m fpr each element
        for j in range(0,3):
            # transform to global
            val_global = dot(transpose(element.T[0:2,0:2]), array((linspace(0, element.length, num=51), val_local[j])))
            x_global = val_global[0]
            val_global = val_global[1]          

            x_def = zeros(51)
            y_def = zeros(51)

            for i in range(0,51):
                x_def[i] = x_1 + x_global[i]
                y_def[i] = y_1 + val_global[i]
        
            sub[j].plot(x_org, y_org,'k')
            sub[j].plot(x_def, y_def,'r')            

        #plt.fill_between(x_def, y_org, y_def)
        x.extend([x_1, x_2])   #x and y vector with all coordinates
        y.extend([y_1, y_2])
      
    #set Limits
    xmin = math.floor(min(x))   # set min and max values
    xmax = math.ceil(max(x))
    xmin = xmin - math.ceil((xmax-xmin)/10)
    xmax = xmax + math.ceil((xmax-xmin)/10)
    ymin = math.floor(min(y))
    ymax = math.ceil(max(y))
    ymin = ymin - math.ceil((ymax-ymin)/10)
    ymax = ymax + math.ceil((ymax-ymin)/10)

    #set settings for each subplot
    for j in range(0,3):
        sub[j].set_xlim( xmin, xmax )   # set axes
        sub[j].set_ylim( ymax, ymin )
        sub[j].set_xlabel('X axis')
        sub[j].set_ylabel('Y axis')

        sub[j].set_xticks(np.arange(xmin, xmax, 1)) # set grid
        sub[j].set_yticks(np.arange(ymin, ymax, 1))
        sub[j].grid()

        sub[j].scatter(x, y, marker='o', c='r', s=5)   # draw nodes
        sub[j].axis('equal')


    #plt.savefig(name)
    plt.show()
    plt.close()



