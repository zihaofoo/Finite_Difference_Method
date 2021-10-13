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
import time

# def ChannelFlow(N_xi, N_eta, bb, hh, ll):

N_xi = 21
N_eta = 21 
bb = 0.5        # in x domain
hh = 1.0        # in y domain
ll = 3.0
aa = np.sqrt(0.25 * ((ll - bb)**2) - (hh**2))
# time1 = time.time()

# Inputs:
# N_xi : Number Of nodes in the xi direction.
# N_eta: Number of nodes int he eta directionpython3 --version
# bb : The length of the base
# hh : the height of the channel

# Outputs:
# Solution : A xi-eta matrix containing the solution
# Flowrate : The flowrate throught eh channel
# I_xx     : The moment of inertia of the channel

d_xi  = 1.0/(N_xi-1);            # delta xi                     
d_eta = 1.0/(N_eta-1);           # delta eta                      
NumNodes = N_xi * N_eta;   

# Initializing sparse matrices 'A' and 'Node' (coefficient value at i,j)
A = sparse.lil_matrix((NumNodes, NumNodes), dtype=float)
Node = np.arange(0, NumNodes, 1, dtype=float)
Node = Node.reshape((N_xi, N_eta)).T
RHS = np.zeros((NumNodes,1), dtype=float)       # Vector of b (nx x ny, 1)
# print(Node)

# Constructing the Jacobian at each i,j
J = np.zeros((N_xi, N_eta), dtype=float)
for i1 in range(N_xi): 
    for i2 in range(N_eta):  
        J[i1,i2] = hh * (0.5 * bb + aa * i2 * d_eta)

# Interior points of matrix A 
for i1 in range(1,N_xi-1):
    for i2 in range(1,N_eta-1):
        # Point (xi, eta) = (i1, i2)
        # Iterate through all interior points
        ANode_i = Node[i1,i2];                       # Setting A_Matrix position for node i,j   
        ANode_i = int(ANode_i)

        xi_local = i1 * d_xi
        eta_local = i2 * d_eta

        # Local coefficients @ (i1, i2)
        c1 = aa**2.0 * xi_local**2.0 + (hh**2.0)
        c2 = (0.5 * aa * bb * xi_local) + ((aa**2) * eta_local * xi_local)
        c3 = ((0.5 * bb) + (aa * eta_local))**2
        c4 = 0
        c5 = 2.0 * (aa**2) * xi_local
 
        # Insert coefficients for the finite difference approximation
        A[ANode_i, int(Node[i1-1, i2-1])] = -0.5 * c2 / (d_eta * d_xi)
        A[ANode_i, int(Node[i1-1, i2])] = (c1 / (d_eta * d_eta)) - (0.5 * c5 / d_xi)
        A[ANode_i, int(Node[i1-1, i2+1])] = 0.5 * c2 / (d_eta * d_xi)
        A[ANode_i, int(Node[i1, i2-1])] = c3 / (d_eta * d_eta)
        A[ANode_i, int(Node[i1, i2])] = (-2.0 * c1 / (d_xi * d_xi)) + (-2.0 * c3 / (d_eta * d_eta))
        A[ANode_i, int(Node[i1, i2+1])] =  c3 / (d_xi * d_xi)
        A[ANode_i, int(Node[i1+1, i2-1])] = 0.5 * c2 / (d_eta * d_xi)
        A[ANode_i, int(Node[i1+1, i2])] =  (c1 / (d_eta * d_eta)) + (0.5 * c5 / d_xi)
        A[ANode_i, int(Node[i1+1, i2+1])] = -0.5 * c2 / (d_eta * d_xi)
        RHS[ANode_i] = -(J[i1,i2])**2.0     # Equation on RHS

# Bottom domain B'C'
for j1 in range(1,N_xi-1):
    j2 = 0                           # First row at the bottom
    ANode_i = int(Node[j1,j2])       # Setting A_Matrix position for node i,j   
    A[ANode_i, int(Node[j1,j2])] = 1.0

# Top domain A'D'
for j1 in range(1,N_xi-1):    
    j2 = N_eta - 1
    ANode_i = int(Node[j1,j2])       # Setting A_Matrix position for node i,j   

    if j1 == N_xi - 2:
        A[ANode_i, int(Node[j1-1,j2])] = (0.5 * (j1 - 1) * aa) 
        A[ANode_i, int(Node[j1,j2])] = (0.5 * bb + aa) * (1.5 / d_eta) 
        A[ANode_i, int(Node[j1+1,j2])] = (-0.5 * (j1 - 1) * aa) 
        A[ANode_i, int(Node[j1,j2-1])] = (0.5 * bb + aa) * (-2.0 / d_eta) 
        A[ANode_i, int(Node[j1,j2-2])] = (0.5 * bb + aa) * (0.5 / d_eta)
    else:
        A[ANode_i, int(Node[j1,j2])] = (1.5 * (j1 - 1) * aa ) + ((0.5 * bb + aa) * (1.5 / d_eta)) 
        A[ANode_i, int(Node[j1+1,j2])] = -2.0 * aa * (j1-1) 
        A[ANode_i, int(Node[j1+2,j2])] = 0.5 * (j1-1) * aa 
        A[ANode_i, int(Node[j1,j2-1])] = (0.5 * bb + aa) * (-2.0 / d_eta) 
        A[ANode_i, int(Node[j1,j2-2])] = (0.5 * bb + aa) * (0.5 / d_eta) 

# Left domain A'B'
for j2 in range(1,N_eta-1):    
    j1 = 0
    ANode_i = int(Node[j1,j2])       # Setting A_Matrix position for node i,j   
    A[ANode_i, int(Node[j1,j2])] = -1.5 / d_xi 
    A[ANode_i, int(Node[j1+1,j2])] = 2.0 / d_xi 
    A[ANode_i, int(Node[j1+2,j2])] = -0.5 / d_xi 

# Right domain C'D'
for j2 in range(1,N_eta-1):
    j1 = N_xi - 1
    ANode_i = int(Node[j1,j2])       # Setting A_Matrix position for node i,j   
    A[ANode_i, int(Node[j1,j2])] = 1.0

# Bottom left corner point
ANode_i = int(Node[0,0])
A[ANode_i, int(Node[0,0])] = 1.0

# Bottom right corner point
ANode_i = int(Node[N_xi-1,0])
A[ANode_i, int(Node[N_xi-1,0])] = 1.0

# Top left corner point
lfcorner_x = 0
lfcorner_y = N_eta-1
ANode_i = int(Node[lfcorner_x, lfcorner_y ])
A[ANode_i, int(Node[lfcorner_x, lfcorner_y ])] = (1.5 * (lfcorner_x - 1) * aa ) + ((0.5 * bb + aa) * (1.5 / d_eta))
A[ANode_i, int(Node[lfcorner_x+1, lfcorner_y ])] = -2.0 * aa * (lfcorner_x-1) 
A[ANode_i, int(Node[lfcorner_x+2, lfcorner_y ])] = 0.5 * (lfcorner_x-1) * aa 
A[ANode_i, int(Node[lfcorner_x, lfcorner_y-1])] = (0.5 * bb + aa) * (-2.0 / d_eta) 
A[ANode_i, int(Node[lfcorner_x, lfcorner_y-2])] = (0.5 * bb + aa) * (0.5 / d_eta) 

# Top right corner point
ANode_i = int(Node[N_xi-1,N_eta-1])
A[ANode_i, int(Node[N_xi-1,N_eta-1])] = 1.0
# np.savetxt("MatrixA.csv", A, delimiter=",")
# print(A, RHS) 

u_sol = linalg.spsolve(A,RHS)
u_sol = u_sol.reshape((N_xi, N_eta))
# print("U_sol:", u_sol)
Q = 0.0

for k1 in range(N_xi):
    for k2 in range(N_eta):
        Q = Q + d_xi * d_eta * u_sol[k1, k2] * hh * (0.5 * bb + aa * d_eta * (k2-1))
print("Q is:", Q)

# time2 = time.time()
# print("Duration:", time2-time1)

xi_vec = np.arange(0, N_xi, 1) * d_xi
eta_vec = np.arange(0, N_xi, 1) * d_xi

xi_grid, xi_grid_T = np.meshgrid(xi_vec, xi_vec)
eta_grid_T, eta_grid = np.meshgrid(eta_vec, eta_vec)

# print(xi_grid, eta_grid)

x_grid = 0.5 * bb * xi_grid + aa * xi_grid * eta_grid
y_grid = hh * eta_grid

# print(x_grid, y_grid)
color_ticks = np.arange(0, 0.22, 0.02)
fig, ax = plt.subplots()
cs = ax.contour(x_grid, y_grid, u_sol, levels=color_ticks, colors='red', linestyles='dashed')
plt.clabel(cs, inline=True, fontsize=10)
cp = ax.contourf(x_grid, y_grid, u_sol, levels=color_ticks)
cbar = fig.colorbar(cp)
plt.savefig("Fig1D.png", dpi=600)
plt.show()

# return Solution, FlowRate, I_xx

