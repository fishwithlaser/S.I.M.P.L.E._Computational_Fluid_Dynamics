############## 
# title:  Using Upwind and Central Difference to solve linear and vortex flow
# by: thomas kosciuch
# email: thomas.kosciuch@mail.utoronto.ca
#
# note, plots save in a 'plot' directory, 
# and text files save in a 'txt' directory
#
# SOURCES:
# Math: https://github.com/mitya57/python-markdown-math
# Scipy: https://github.com/scipy/scipy
# Time: https://github.com/time
# Sys: https://github.com/sys
# tqdm: https://github.com/tqdm/ 
# numpy: https://github.com/numpy/numpy
# future: https://github.com/python/cpython/blob/2.7/Lib/__future__.py
# matplotlib: https://github.com/matplotlib/matplotlib
#
# __ I M P O R T I N G #
from __future__ import division    # solution to complete division using floats
import numpy as np                 # 
import numpy.linalg                # direct solver for solution
import math                        # 
import scipy.sparse.linalg         # incase we want our nodes to be > 160

import time                        # enables timing for debugging purposes
import sys                         # allows integer inputs
from tqdm import tqdm              # enables asking questions
import matrixPlotter as Plotter
import buildMatrix as matrix
#  __  K N O W N S  _  #
# note node numbers, u's and gammas all subject to change
x_length  = 1                          # length of X, float 
y_length  = 1                          # length of Y, float
x_nodes   = 40                         # number of nodes, integer
y_nodes   = 40                         # number of nodes, integer
Circ      = 'y'
Scheme    = 'u'
debug     = 'n' 
x_u       = 0.1                        #
x_v       = 0.1                        #
Gamma     = 0.1                        # diffusion coefficient
phi_N     = 100                        # boundary at N 
phi_E     = 0                          # boundary at E
phi_S     = 0                          # boundary at N 
phi_W     = 100                        # boundary at E
density   = 1                          # rho term
Boundary = [phi_W, phi_E, phi_S, phi_N]

# SIMPLE MATH
dim = y_nodes *  x_nodes               # dimension of a-matrix
dx  = x_length / x_nodes               # node size for static grid along x
dy  = y_length / y_nodes               # node size for static gris along y
x_pos, y_pos = matrix.meshposition(x_nodes,y_nodes,dx,dy)    

#   C A L C U L A T I N G   V E L O C I T Y   #

x_u   = [x_u] * dim                    # u(x) values are populated
x_v   = [x_v] * dim                    # u(y) values
                                       # vortex will over-ride:
if Circ == "y":
    print("calculating circular velocities")
    matrix.CircV(x_pos,y_pos,x_length,y_length, dim)

# C A L C U L A T I N G    P E C L E T   No. #
F_x     = [0] * dim
F_y     = [0] * dim
D_x = Gamma / dx
D_y = Gamma / dy                    # This allows the calculation of a 1-D system

if x_nodes == 1:                    # along x- and y- 
    D_x = 0                         # which is very helpful for checking programming
if y_nodes == 1:
    D_y = 0

for i in range(dim):
    F_x[i] = density * x_u[i]       # calculates F-values along x
    F_y[i] = density * x_v[i]       # calculates F-values along y
   
# Central Difference
if Scheme == "c":
    a, B = matrix.CDiff(x_nodes, y_nodes, D_x, D_y, F_x, F_y, Boundary, debug)
    print("central difference")

# Upwind Difference
if Scheme == "u":
    a, B = matrix.UDiff(x_nodes, y_nodes, D_x, D_y, F_x, F_y, Boundary, debug)
    print("upwind difference")


xpos = [0] * x_nodes                      # calculating position along x
ypos = [0] * y_nodes                      # calculating position along y
xpos[0] = dx / 2                          # populating the first x-node
ypos[0] = dy / 2                          # populating the first y-node

for i in range(x_nodes-1):                
    xpos[i+1] = xpos[i] + dx
for i in range(y_nodes-1):
    ypos[i+1] = ypos[i] + dy

if x_nodes < 170:                        # solution is direct if fewer than 170 nodes
    print("solving - direct")
    X = np.linalg.solve(a, B)
    print("direct - successful")
else:                                    # iterative solution for large matricies
    print("solving - Iterative")
    a_temp = np.asmatrix(a)
    x_temp = scipy.sparse.linalg.gmres(a_temp,B)
    if x_temp[1] is 0:
        print("interative - success") 
    else:
        print("iterative - failure, I'm sorry - goodbye")
        quit
    X = x_temp[0]


Ans = matrix.List2Matrix(X,x_nodes,y_nodes)
print(Ans)

# creating name for plotter
name = "dif=" + str(Gamma) + " Sch=" + Scheme + " N=" + str(dim) + " u1:" + str(round(x_u[1],2))

# graphing plot and error plot
Plotter.Contour(Ans, B, name=name)
Plotter.Error(Ans, name=name)


