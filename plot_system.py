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


def PlotSystem(Nodes, Elements, variable, name, scale):
    nnodes = len(Nodes)
    nelements = len(Elements)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    x = []  # preallocation
    y = []
    x_val = []  # preallocation
    y_val = []

    for element in Elements:
        # original coordinates -> 2D Line inbetween
        x_1 = element.geometry.nodes[0].coordinates[0]
        y_1 = element.geometry.nodes[0].coordinates[1]
        x_2 = element.geometry.nodes[1].coordinates[0]
        y_2 = element.geometry.nodes[1].coordinates[1]
        line = Line2D([x_1,x_2],[y_1,y_2])
        line.set_linestyle('dashed')
        line.set_color('gray')
        ax.add_line(line)   # add line2D

        # deflection curves of elements
        [x_def, y_def] = element.CalculateDeflectionCurve(scale)
        plt.plot(x_def,y_def,'k')

        x.extend([x_1, x_2])   #x and y vector with all coordinates
        y.extend([y_1, y_2])
        x_val.extend(x_def)   #val vector für min and max
        y_val.extend(y_def)

    xmin = math.floor(min(x_val))   # set min and max values
    xmax = math.ceil(max(x_val))
    xmin = xmin - math.ceil((xmax-xmin)/10)
    xmax = xmax + math.ceil((xmax-xmin)/10)
    ymin = math.floor(min(y_val))
    ymax = math.ceil(max(y_val))
    ymin = ymin - math.ceil((ymax-ymin)/10)
    ymax = ymax + math.ceil((ymax-ymin)/10)


    ax.set_xlim( xmin, xmax )   # set axes
    ax.set_ylim( ymax, ymin )
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')

    ax.set_xticks(np.arange(xmin, xmax, 1)) # set grid
    ax.set_yticks(np.arange(ymin, ymax, 1))
    plt.grid()

    plt.scatter(x, y, marker='o', c='r', s=5)   # draw nodes
    plt.title('%d nodes, %d elements, scale = %d' % (nnodes, nelements, scale) )    # set title
    plt.axes().set_aspect('equal')  # equalize axes

    plt.savefig(name)
    plt.show()
    plt.close()