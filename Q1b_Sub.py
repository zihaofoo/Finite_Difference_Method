# Massachusetts Institute of Technology
# Department of Mechanical Engineering
author = "Zi Hao Foo"
date = "Oct 2021"

import numpy as np
import scipy as sp
from scipy.sparse import linalg
from scipy import sparse
import matplotlib.pyplot as plt
import pandas as pd

N_x = 7
N_y = 7 
d_x  = 1.0/(N_x-1);            # delta x                     
d_y = 1.0/(N_y-1);           # delta y                     
NumNodes = N_x * N_y;   

NumNodes_unitlblock = int((N_x - 1) / 6)

# Initializing sparse matrices 'A' and 'Node' (coefficient value at i,j)
A = sparse.lil_matrix((NumNodes, NumNodes), dtype=float)
Node = np.arange(0, NumNodes, 1, dtype=float)
Node = Node.reshape((N_x, N_y)).T
RHS = np.zeros((NumNodes,1), dtype=float)       # Vector of f (nx x ny, 1)

# Interior points of matrix A 
    for i1 in range(1,N_x-1):
        for i2 in range(1,N_y-1):
            # Point (x, y) = (i1, i2)
            # Iterate through all interior points
            ANode_i = Node[i1,i2];                       # Setting A_Matrix position for node i,j   
            ANode_i = int(ANode_i)

            x_local = i1 * d_x
            y_local = i2 * d_y
    
            # Insert coefficients for the finite difference approximation
            A[ANode_i, int(Node[i1-1, i2-1])] = 
            A[ANode_i, int(Node[i1-1, i2])] = 
            A[ANode_i, int(Node[i1-1, i2+1])] = 
            A[ANode_i, int(Node[i1, i2-1])] = 
            A[ANode_i, int(Node[i1, i2])] = 
            A[ANode_i, int(Node[i1, i2+1])] =  
            A[ANode_i, int(Node[i1+1, i2-1])] = 
            A[ANode_i, int(Node[i1+1, i2])] =  
            A[ANode_i, int(Node[i1+1, i2+1])] = 
            RHS[ANode_i] = -(J[i1,i2])**2.0     # Equation on RHS