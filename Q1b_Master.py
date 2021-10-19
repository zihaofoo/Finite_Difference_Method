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
from Q1b_Sub import Solver_Iter, Flux, I_down, Plotting

Source_block = np.array([1, 7, 14, 16])
# N_x = np.array([7, 13, 19, 25])                    # Number of nodes
# N_y = np.array([7, 13, 19, 25])                   # Number of nodes
N_Iter = 5000
N_sol = np.arange(N_Iter)
N_x = np.array([25], dtype=int)                    # Number of nodes
N_y = np.array([25], dtype=int)                   # Number of nodes
# omega = np.arange(1, 1.05, 0.05) 
omega = np.array([1.0])
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
u_sol, err_sol = Solver_Iter(Source_block, int(N_x[0]), int(N_y[0]), omega[0], method1, N_Iter) 
# flux_top, flux_bottom, flux_left, flux_right = Flux(u_sol, int(N_x[0]), int(N_y[0]))
# flux_net = np.concatenate((flux_left, flux_right, flux_bottom, flux_top), axis=0)
# np.savetxt("flux.csv", flux_net)

u_sol_coarse = I_down(u_sol, N_x, N_y)

Plotting(u_sol, N_x, N_y)
Plotting(u_sol_coarse, 13, 13)

np.savetxt('u_sol_fine.csv', u_sol)
np.savetxt('u_sol_coarse.csv', u_sol_coarse)

