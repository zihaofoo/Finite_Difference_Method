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
N_Iter = 100
N_sol = np.arange(N_Iter)
N_x = int(49)                    # Number of nodes
N_y = int(49)               # Number of nodes
# omega = np.arange(1, 1.05, 0.05) 
omega = np.array([0.5, 0.8])
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

## Testing prolongation and restriction matrix
"""
u_sol_fine = I_up(u_sol, N_x, N_y)

Plotting(u_sol, N_x, N_y)
Plotting(u_sol_fine, 49, 49)

np.savetxt('u_sol_1.csv', u_sol)
np.savetxt('u_sol_2.csv', u_sol_fine)
"""

## Multigrid solver
nu_1 = int(1)
nu_2 = int(1)
nu_c = int(4)
# N_factor = int(1)
N_factor = np.array([4], dtype=int)

# Error plot
fig, ax = plt.subplots()

"""
for i1 in range(omega.shape[0]):
    u_sol, err_sol = MultiGrid(Source_block, N_x, N_y, omega[i1], method1, nu_1, nu_2, nu_c, N_factor, N_Iter)
    cp = ax.plot(N_sol, err_sol, label = ('Omega = ' + str(omega[i1])))
"""

for i2 in range(N_factor.shape[0]):
    u_sol, err_sol = MultiGrid(Source_block, N_x, N_y, omega[1], method1, nu_1, nu_2, nu_c, int(N_factor[i2]), N_Iter)
    print(i2)
    cp = ax.plot(N_sol, err_sol, label = ('Omega = ' + str(omega[1])))

ax.set_yscale('log')
ax.set_title('Method = ' + method1)
# ax.set_ybound(lower = 1E-12, upper = 1E+08)
plt.legend()
# plt.savefig(fname = (method1 + '_Iter.png'), dpi=600)
plt.show()
