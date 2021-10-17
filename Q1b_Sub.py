# Massachusetts Institute of Technology
# Department of Mechanical Engineering
author = "Zi Hao Foo"
date = "Oct 2021"

import numpy as np
import scipy as sp
import scipy.linalg as linalg
import scipy.sparse.linalg as splinalg
from scipy import sparse
import matplotlib.pyplot as plt
import pandas as pd

def Solver_Iter(Source_block, N_x, N_y, omega, method, N_iter): 
    # Source_block = np.array([1, 6, 14, 16])
    # N_x = 25                    # Number of nodes
    # N_y = 25                    # Number of nodes
    # omega = 0.5 
    # method = 'Jacobi'

    d_x  = 1.0/(N_x-1);         # delta x                     
    d_y = 1.0/(N_y-1);          # delta y                     
    NumNodes = N_x * N_y;   
    NumNodes_interior = int((N_x - 1) / 6) - 1    # Number of interior nodes per block 

    # Initializing sparse matrices 'A' and 'Node' (coefficient value at i,j)
    A = sparse.lil_matrix((NumNodes, NumNodes), dtype=float)
    Node = np.arange(0, NumNodes, 1, dtype=float)
    Node = Node.reshape((N_x, N_y)).T

    # print(Node[2,0])
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
            A[ANode_i, int(Node[i1-1, i2])] = - 1.0 / (d_x**2.0)
            A[ANode_i, int(Node[i1, i2-1])] = - 1.0 / (d_y**2.0)
            A[ANode_i, int(Node[i1, i2])] = (2.0 / (d_x**2.0)) + (2.0 / (d_y**2.0))
            A[ANode_i, int(Node[i1, i2+1])] = - 1.0 / (d_y**2.0)
            A[ANode_i, int(Node[i1+1, i2])] = - 1.0 / (d_x**2.0) 

    # Bottom domain B'C'
    for j1 in range(1,N_x):
        j2 = 0                           # First row at the bottom
        ANode_i = int(Node[j1,j2])       # Setting A_Matrix position for node i,j   
        A[ANode_i, int(Node[j1,j2])] = 1/(d_x**2.0)# 1.0 

    # Top domain A'D'
    for j1 in range(1,N_x-1):    
        j2 = N_y - 1
        ANode_i = int(Node[j1,j2])       # Setting A_Matrix position for node i,j  
        A[ANode_i, int(Node[j1,j2])] = 1/(d_x**2.0) # 1.0 

    # Left domain A'B'
    for j2 in range(1,N_y-1):    
        j1 = 0
        ANode_i = int(Node[j1,j2])       # Setting A_Matrix position for node i,j  
        A[ANode_i, int(Node[j1,j2])] = 1/(d_x**2.0) # 1.0 

    # Right domain C'D'
    for j2 in range(1,N_y-1):
        j1 = N_x - 1
        ANode_i = int(Node[j1,j2])       # Setting A_Matrix position for node i,j   
        A[ANode_i, int(Node[j1,j2])] = 1/(d_x**2.0) # 1.0

    # Bottom left corner point
    ANode_i = int(Node[0,0])
    A[ANode_i, int(Node[0,0])] = 1/(d_x**2.0) # 1.0

    # Bottom right corner point
    ANode_i = int(Node[N_x-1,0])
    A[ANode_i, int(Node[N_x-1,0])] = 1/(d_x**2.0) # 1.0

    # Top left corner point
    ANode_i = int(Node[0, N_y-1])
    A[ANode_i, int(Node[0,N_y-1])] = 1/(d_x**2.0) # 1.0

    # Top right corner point
    ANode_i = int(Node[N_x-1,N_y-1])
    A[ANode_i, int(Node[N_x-1,N_y-1])] = 1/(d_x**2.0) # 1.0

    # Vector RHS 
    column_block = np.mod(Source_block - 1, 4) + 1          # Column index of source blocks
    row_block = np.floor_divide(Source_block - 1, 4) + 1    # Row index of source blocks

    for s1 in range(Source_block.shape[0]):
        btm_left_x = row_block[s1] * (NumNodes_interior + 1)
        btm_left_y = column_block[s1] * (NumNodes_interior + 1)
        top_right_x = btm_left_x + NumNodes_interior + 1
        top_right_y = btm_left_y + NumNodes_interior + 1

        for i1 in range(btm_left_x, top_right_x+1):
            for i2 in range(btm_left_y, top_right_y+1):
                # Point (x, y) = (i1, i2)
                # Iterate through all interior points
                ANode_i = Node[i1,i2];                       # Setting A_Matrix position for node i,j   
                ANode_i = int(ANode_i)
                # Insert coefficients for the RHS vector
                RHS[ANode_i] = 1.0   # Equation on RHS

    D = np.diag(np.diag(A.toarray()))
    L = np.tril(-A.toarray(), k=-1)
    U = np.triu(-A.toarray(), k=1)
    I = np.identity(NumNodes, dtype=float)

    # N_iter = 10000
    u_sol = np.zeros((NumNodes,1), dtype=float)       # Vector of f (nx x ny, 1)

    # omega = 2.0/3
    # method = 'Jacobi'
    # method = 'Gauss'

    err_sol = np.zeros(N_iter)
    if method == 'Jacobi':
        D_inv = linalg.inv(D)
        R = I - (D_inv @ A)    
        f = D_inv @ RHS
    elif method == 'Gauss':
        D_L_inv = linalg.inv(D - L)
        R = D_L_inv @ U   
        f = D_L_inv @ RHS

    for t1 in range(N_iter):
        u_sol = (omega * ((R @ u_sol) + f)) + ((1.0 - omega) * u_sol)
        err_sol[t1] = linalg.norm((A @ u_sol - RHS), ord=2) 

    # Linear sparse solver (for checking only)
    # u_sol = splinalg.spsolve(A,RHS)

    ## Plotting
    """
    u_sol_matrix = u_sol.reshape((N_x, N_y))
    x_vec = np.arange(0, N_x, 1) * d_x
    y_vec = np.arange(0, N_y, 1) * d_y
    fig, ax = plt.subplots()
    cp = ax.contourf(x_vec, y_vec, u_sol_matrix)
    cbar = fig.colorbar(cp)
    plt.show()
    """
    return u_sol, err_sol

