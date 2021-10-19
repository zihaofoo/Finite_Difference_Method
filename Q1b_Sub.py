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
    A = sparse.lil_matrix((NumNodes, NumNodes), dtype=np.double)
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
    u_sol = np.zeros((NumNodes,1), dtype=float)   # Vector of f (nx x ny, 1)

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
        err_sol[t1] = linalg.norm(np.double(A @ u_sol - RHS), ord=2) 

    # Linear sparse solver (for checking only)
    """
    u_sol_linear = splinalg.spsolve(A,RHS)
    u_sol_linear = np.reshape(u_sol_linear,(u_sol_linear.size, 1))
    err_linear = u_sol - u_sol_linear
    print(linalg.norm(err_linear, ord=2))
    """

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

def Plotting(u_sol, N_x, N_y):
    N_x = int(N_x)
    N_y = int(N_y)
    d_x  = 1.0/(N_x-1);         # delta x                     
    d_y = 1.0/(N_y-1);          # delta y  
    ## Plotting
    u_sol_matrix = u_sol.reshape((N_x, N_y))
    x_vec = np.arange(0, N_x, 1) * d_x
    y_vec = np.arange(0, N_y, 1) * d_y
    fig, ax = plt.subplots()
    cp = ax.contourf(x_vec, y_vec, u_sol_matrix)
    cbar = fig.colorbar(cp)
    plt.show()

    return 

'''
function uh = MG_rec( uh, fh, nu1, nu2, omega, Nmin)
    uh = relax( uh, fh, nu1, omega);

    if length(uh) > Nmin
        rh = residual(uh, fh);
        r2h = restrict(rh);
        e2h = MG_rec( 0*r2h, r2h, nu1, nu2, omega, Nmin);
        eh = prolongate(e2h);
        uh = uh + eh;
    end
    
    uh = relax( uh, fh, nu2, omega);
end
'''

def MultiGrid(Source_block, N_x, N_y, omega, method, nu_1, nu_2, N_factor, N_iter):
    N_temp = N_x
    N_x = int(N_x)
    N_y = int(N_y)

    for i1 in range(N_factor):
        N_temp = (N_temp + 1) / 2
        N_min = int(N_temp)
    
    NumNodes = N_x * N_y;   
    u_sol = np.zeros((NumNodes,1), dtype=float)   # Vector of f (nx x ny, 1)
    err_sol = np.zeros(N_iter, dtype=float)
    
    for t1 in range(N_iter):
        u_sol = MultiGrid_Sub(u_sol, Source_block, N_x, N_y, omega, method, nu_1, nu_2, N_min)
        res_sol = Solver_Relax(u_sol, Source_block, N_x, N_y, omega, method, N_iter=1, mode='Residual')
        err_sol[t1] = linalg.norm(res_sol, ord=2) 
        # print(t1)

    return u_sol, err_sol

def MultiGrid_Sub(u_sol, Source_block, N_x, N_y, omega, method, nu_1, nu_2, N_min): 
    # print('Inception')
    u_sol = Solver_Relax(u_sol, Source_block, N_x, N_y, omega, method, N_iter=nu_1)      # Relax by nu 1 times

    if N_x > N_min: 
        res_h = Solver_Relax(u_sol, Source_block, N_x, N_y, omega, method, N_iter=1, mode='Residual')   # Calculate residual
        print(linalg.norm(res_h, ord=2))
        N_x_coarse = int((N_x + 1) / 2)     # N_x of coarse grid
        N_y_coarse = int((N_y + 1) / 2)     # N_y of coarse grid
        res_2h = I_down(res_h, N_x, N_y, skip_factor=int(2)) 
        # print("Pre-Inception")
        # print(res_h.shape)
        # print(res_2h.shape, N_x_coarse)
        err_2h = MultiGrid_Sub(0 * res_2h, Source_block, N_x_coarse, N_y_coarse, omega, method, nu_1, nu_2, N_min)
        err_h = I_up(err_2h, N_x_coarse, N_y_coarse, skip_factor=int(2))
        u_sol = u_sol + err_h

    u_sol = Solver_Relax(u_sol, Source_block, N_x, N_y, omega, method, N_iter=nu_2)      # Relax by nu 2 times

    return u_sol

def I_up(u_sol, N_x, N_y, skip_factor=int(2)): 
    # N_x, N_y are dimensions of u_sol (coarse)
    # u_sol on skip_factor * h, u_sol_fine on h
    N_x = int(N_x)        # dimension of coarse mesh
    N_y = int(N_y)        # dimension of coarse mesh                  
    NumNodes = N_x * N_y;   
    Node = np.arange(0, NumNodes, 1, dtype=float)
    Node = Node.reshape((N_x, N_y)).T

    u_sol_fine = []

    for k1 in range(0, skip_factor * N_x - 1, 1):
        for k2 in range(0, skip_factor * N_y - 1, 1):
            if np.mod(k1, 2) == 0:
                if np.mod(k2, 2) == 0:
                    ANode_i = Node[int(k2/2), int(k1/2)];                       # Setting A_Matrix position for node i,j   
                    ANode_i = int(ANode_i)
                    u_sol_fine.append(u_sol[ANode_i])  # h, {2i, 2j} = 2h, {i,j}
                else: 
                    u_sol_fine.append((0.5 * u_sol[int(Node[ int((k2-1)/2), int(k1/2)] )]) + (0.5 * u_sol[int(Node[ int((k2+1)/2), int(k1/2)] )])) # h, {2i, 2j+1} 
            else:
                if np.mod(k2, 2) == 0:
                    u_sol_fine.append((0.5 * u_sol[int(Node[ int(k2/2), int((k1-1)/2)]) ]) + (0.5 * u_sol[int(Node[ int(k2/2), int((k1+1)/2) ])])) # h, {2i+1, 2j} 
                else:
                    u_sol_fine.append((0.25 * u_sol[int(Node[ int(k2/2), int(k1/2)]) ]) + (0.25 * u_sol[int(Node[ int(k2/2), int((k1+1)/2)]) ]) + (0.25 * u_sol[int(Node[ int((k2+1)/2), int(k1/2)]) ]) + (0.25 * u_sol[int(Node[ int((k2+1)/2), int((k1+1)/2)])] )) # h, {2i+1, 2j+1}

    u_sol_fine = np.array(u_sol_fine)
    
    return u_sol_fine

def I_down(u_sol, N_x, N_y, skip_factor=int(2)):
    # N_x, N_y are dimensions of u_sol (fine)
    # u_sol on h, u_sol_coarse on skip_factor * h
    N_x = int(N_x)
    N_y = int(N_y)                
    NumNodes = N_x * N_y;   
    Node = np.arange(0, NumNodes, 1, dtype=float)
    Node = Node.reshape((N_x, N_y)).T

    u_sol_coarse = []

    for k1 in range(0, N_x, skip_factor):
        for k2 in range(0, N_y, skip_factor):
            ANode_i = Node[k2,k1];                       # Setting A_Matrix position for node i,j   
            ANode_i = int(ANode_i)
            u_sol_coarse.append(u_sol[ANode_i])
    u_sol_coarse = np.array(u_sol_coarse)

    return u_sol_coarse

def Solver_Relax(u_sol_initial, Source_block, N_x, N_y, omega, method, N_iter, mode='Relax'): 
    d_x  = 1.0/(N_x-1);         # delta x                     
    d_y = 1.0/(N_y-1);          # delta y                     
    NumNodes = N_x * N_y;   
    NumNodes_interior = int((N_x - 1) / 6) - 1    # Number of interior nodes per block 

    # Initializing sparse matrices 'A' and 'Node' (coefficient value at i,j)
    A = sparse.lil_matrix((NumNodes, NumNodes), dtype=np.double)
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

    u_sol = u_sol_initial  # Vector of f (nx x ny, 1)

    if mode == "Residual":
        return (RHS - A @ u_sol)

    D = np.diag(np.diag(A.toarray()))
    L = np.tril(-A.toarray(), k=-1)
    U = np.triu(-A.toarray(), k=1)
    I = np.identity(NumNodes, dtype=float)

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

    return u_sol


def Flux(u_sol, N_x, N_y): 
    d_x  = 1.0/(N_x-1);         # delta x                     
    d_y = 1.0/(N_y-1);          # delta y      
    NumNodes = N_x * N_y;   
    Node = np.arange(0, NumNodes, 1, dtype=float)
    Node = Node.reshape((N_x, N_y)).T

    # Bottom domain B'C'
    flux_bottom = np.zeros(N_x, dtype=float)
    for j1 in range(1,N_x-1):
        j2 = 0                              # First row at the bottom
        ANode_i1 = int(Node[j1,j2])         # Setting A_Matrix position for node i,j   
        ANode_i2 = int(Node[j1,j2+1])       # Setting A_Matrix position for node i,j   
        ANode_i3 = int(Node[j1,j2+2])       # Setting A_Matrix position for node i,j   
        flux_bottom[j1] = -(1.5 / d_y) * u_sol[ANode_i1] + (2.0 / d_y) * u_sol[ANode_i2] - (0.5 / d_y) * u_sol[ANode_i3] 
    # print(j1, "Bottom")

    flux_top = np.zeros(N_x, dtype=float)
    # Top domain A'D'
    for j1 in range(1,N_x-1):    
        j2 = N_y - 1
        ANode_i1 = int(Node[j1,j2])       # Setting A_Matrix position for node i,j   
        ANode_i2 = int(Node[j1,j2-1])       # Setting A_Matrix position for node i,j   
        ANode_i3 = int(Node[j1,j2-2])       # Setting A_Matrix position for node i,j   
        flux_top[j1] = (1.5 / d_y) * u_sol[ANode_i1] - (2.0 / d_y) * u_sol[ANode_i2] + (0.5 / d_y) * u_sol[ANode_i3] 
    # print(j1, "Top")

    flux_left = np.zeros(N_y, dtype=float)
    # Left domain A'B'
    for j2 in range(1,N_y-1):    
        j1 = 0
        ANode_i1 = int(Node[j1,j2])       # Setting A_Matrix position for node i,j   
        ANode_i2 = int(Node[j1+1,j2])       # Setting A_Matrix position for node i,j   
        ANode_i3 = int(Node[j1+2,j2])       # Setting A_Matrix position for node i,j        
        flux_left[j2] = -(1.5 / d_x) * u_sol[ANode_i1] + (2.0 / d_x) * u_sol[ANode_i2] - (0.5 / d_x) * u_sol[ANode_i3] 
    # print(j2, "Left")

    flux_right = np.zeros(N_y, dtype=float)
    # Right domain C'D'
    for j2 in range(1,N_y-1):
        j1 = N_x - 1
        ANode_i1 = int(Node[j1,j2])       # Setting A_Matrix position for node i,j   
        ANode_i2 = int(Node[j1-1,j2])       # Setting A_Matrix position for node i,j   
        ANode_i3 = int(Node[j1-2,j2])       # Setting A_Matrix position for node i,j    
        flux_right[j2] = (1.5 / d_x) * u_sol[ANode_i1] - (2.0 / d_x) * u_sol[ANode_i2] + (0.5 / d_x) * u_sol[ANode_i3] 
    # print(j2, "Right")

    return flux_top, flux_bottom, flux_left, flux_right