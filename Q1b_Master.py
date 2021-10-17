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
from Q1b_Sub import Solver_Iter

Source_block = np.array([1, 6, 14, 16])
# N_x = np.array([7, 13, 19, 25])                    # Number of nodes
# N_y = np.array([7, 13, 19, 25])                   # Number of nodes
N_Iter = 2000
N_sol = np.arange(N_Iter)
N_x = np.array([25])                    # Number of nodes
N_y = np.array([25])                   # Number of nodes
omega = np.arange(1, 1.9, 0.05) 
method1 = 'Gauss'

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