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
from Q1a_Sub import ChannelFlow

N = np.array([11, 21, 41])
N_ref = 81
# N = np.arange(11, 161, 10)
bb = 1.0        # in x domain
hh = 1.0        # in y domain
ll = 3.0
Q = np.zeros(len(N), dtype=float)
u_sol = np.zeros((len(N), max(N**2)))

Q[0], u_sol_0 = ChannelFlow(N_eta=N[0], N_xi=N[0], bb=bb, hh=hh, ll=ll)
Q[1], u_sol_1 = ChannelFlow(N_eta=N[1], N_xi=N[1], bb=bb, hh=hh, ll=ll)
Q[2], u_sol_2 = ChannelFlow(N_eta=N[2], N_xi=N[2], bb=bb, hh=hh, ll=ll)
# Q[3], u_sol_3 = ChannelFlow(N_eta=N[3], N_xi=N[3], bb=bb, hh=hh, ll=ll)
Q_ref, u_sol_5 = ChannelFlow(N_eta=N_ref, N_xi=N_ref, bb=bb, hh=hh, ll=ll)
u_sol_1_coarse = np.zeros((min(N)**2, min(N)**2))
u_sol_2_coarse = np.zeros((min(N), min(N)))
u_sol_3_coarse = np.zeros((min(N), min(N)))

"""
print(u_sol_1_coarse.shape[0])
k1 = 0
factor_mesh = (N[1] - 1) / (N[0] - 1)
print(u_sol_1.shape[0])
for j1 in range(u_sol_1.shape[0]):
    if np.mod(j1, factor_mesh) == 0:
        print(k1, j1)
        u_sol_1_coarse[k1] = u_sol_1[j1]
        k1 = k1 + 1
        
print(u_sol_1_coarse, u_sol_0)
"""

Q = np.abs(Q - Q_ref)
alp = np.zeros(len(Q)-1)
for i2 in range(len(Q)-1):
    alp[i2] = np.log(Q[i2+1] / Q[i2]) / np.log(N[i2+1] / N[i2])

y_N = N**(-2.0)
fig, ax = plt.subplots()
ax.loglog(N, Q, color='red', label='Q (Convergence)')
ax.loglog(N, y_N, color='black', label='alpha = 2', linestyle='dashed')
ax.set_xlabel('N')
ax.set_ylabel('Q')
ax.set_title('l = ' + str(ll) + ', b = ' + str(bb) + ', h =' + str(hh))
plt.legend()
plt.show()


