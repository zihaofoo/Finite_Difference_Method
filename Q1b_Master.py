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
from Q1b_Sub import I_up, Solver_Iter, Flux, I_up, I_down, Plotting, MultiGrid

Source_block = np.array([1, 7, 14, 16])
# N_x = np.array([7, 13, 19, 25])                    # Number of nodes
# N_y = np.array([7, 13, 19, 25])                   # Number of nodes
N_Iter = 2
N_sol = np.arange(N_Iter)
N_x = int(25)                    # Number of nodes
N_y = int(25)               # Number of nodes
# omega = np.arange(1, 1.05, 0.05) 
omega = np.array([0.5])
method1 = 'Gauss'

## Convergence plot
"""
fig, ax = plt.subplots()
for i1 in range(omega.shape[0]):
    u_sol, err_sol = Solver_Iter(Source_block, int(N_x[0]), int(N_y[0]), omega[i1], method1, N_Iter) 
    cp = ax.plot(N_sol, err_sol, label = ('Omega =' + str(omega[i1])))

ax.set_yscale('log')
ax.set_title('Method = ' + method1)
# ax.set_ybound(lower = 1E-12, upper = 1E+08)
plt.legend()
plt.savefig(fname = (method1 + '_Iter.png'), dpi=600)
plt.show()
"""

## Calculate flux boundary
# u_sol, err_sol = Solver_Iter(Source_block, int(N_x[0]), int(N_y[0]), omega[0], method1, N_Iter) 
# flux_top, flux_bottom, flux_left, flux_right = Flux(u_sol, int(N_x[0]), int(N_y[0]))
# flux_net = np.concatenate((flux_left, flux_right, flux_bottom, flux_top), axis=0)
# np.savetxt("flux.csv", flux_net)

## Testing prolongation and resitrction matrix
"""
u_sol_fine = I_up(u_sol, N_x, N_y)

Plotting(u_sol, N_x, N_y)
Plotting(u_sol_fine, 49, 49)

np.savetxt('u_sol_1.csv', u_sol)
np.savetxt('u_sol_2.csv', u_sol_fine)
"""

## Multigrid solver
method = 'Gauss'
nu_1 = int(1)
nu_2 = int(1)
N_factor = int(2)

u_sol, err_sol = MultiGrid(Source_block, N_x, N_y, omega, method, nu_1, nu_2, N_factor, N_Iter)

# Error plot
fig, ax = plt.subplots()
cp = ax.plot(N_sol, err_sol, label = ('Omega =' + str(omega)))
ax.set_yscale('log')
ax.set_title('Method = ' + method1)
# ax.set_ybound(lower = 1E-12, upper = 1E+08)
plt.legend()
# plt.savefig(fname = (method1 + '_Iter.png'), dpi=600)
plt.show()
