# Massachusetts Institute of Technology
# Department of Mechanical Engineering
author = "Zi Hao Foo"
date = "Oct 2021"

import numpy as np
import scipy as sp
from scipy import linalg as linalg
from scipy.sparse import linalg as splinalg
from scipy import sparse
import matplotlib.pyplot as plt
import pandas as pd
from Q1a_Sub import ChannelFlow

"""
## Part 6
N = np.array([11, 21, 41])
N_base = N[0]
N_ref = 81
# N = np.arange(11, 161, 10)
bb = 0.5      # in x domain
hh = 1.0        # in y domain
ll = 3.0
Q = np.zeros(len(N), dtype=float)
l2_err = np.zeros(len(N), dtype=float)
linf_err = np.zeros(len(N), dtype=float)
u_sol = np.zeros((len(N), max(N**2)))

Q[0], u_sol_0, u_sol_0_coarse = ChannelFlow(N_eta=N[0], N_xi=N[0], bb=bb, hh=hh, ll=ll, N_base=N_base)
Q[1], u_sol_1, u_sol_1_coarse = ChannelFlow(N_eta=N[1], N_xi=N[1], bb=bb, hh=hh, ll=ll, N_base=N_base)
Q[2], u_sol_2, u_sol_2_coarse  = ChannelFlow(N_eta=N[2], N_xi=N[2], bb=bb, hh=hh, ll=ll, N_base=N_base)
Q_ref, u_sol_3, u_sol_3_coarse  = ChannelFlow(N_eta=N_ref, N_xi=N_ref, bb=bb, hh=hh, ll=ll, N_base=N_base)

l2_err[0] = linalg.norm(u_sol_0_coarse - u_sol_3_coarse, ord=2) 
l2_err[1] = linalg.norm(u_sol_1_coarse - u_sol_3_coarse, ord=2)
l2_err[2] = linalg.norm(u_sol_2_coarse - u_sol_3_coarse, ord=2)

linf_err[0] = linalg.norm(u_sol_0_coarse - u_sol_3_coarse, ord=np.inf)
linf_err[1] = linalg.norm(u_sol_1_coarse - u_sol_3_coarse, ord=np.inf)
linf_err[2] = linalg.norm(u_sol_2_coarse - u_sol_3_coarse, ord=np.inf)

Q = np.abs(Q - Q_ref)
alp = np.zeros(len(Q)-1)
for i2 in range(len(Q)-1):
    alp[i2] = np.log(Q[i2+1] / Q[i2]) / np.log(N[i2+1] / N[i2])

print(alp)
y_N = N**(-2.0)
fig, ax = plt.subplots()
ax.loglog(N, Q, color='red', label='Q (Convergence)')
ax.loglog(N, l2_err, color='blue', label='L2 norm (Convergence)')
ax.loglog(N, linf_err, color='maroon', label='Linf norm (Convergence)')
ax.loglog(N, y_N, color='black', label='alpha = 2', linestyle='dashed')
ax.set_xlabel('N')
ax.set_ylabel('Q')
ax.set_title('l = ' + str(ll) + ', b = ' + str(bb) + ', h =' + str(hh))
ax.set_ybound(lower = 1E-01, upper = 1E-05)
plt.legend()
plt.savefig(fname= str(bb) + 'convergence.PDF')
plt.show()

"""
# Part 7
N = np.array([21])
ll = 3.0
bb = np.arange(0.0, 1.1, 0.1)        # in x domain
hh = np.arange(0.1, 1.1, 0.1)        # in y domain
N_bb, N_hh = (bb.shape[0], hh.shape[0])
Q = np.zeros((N_bb, N_hh), dtype=float)
I_xx = np.zeros((N_bb, N_hh), dtype=float)

for i1 in range(N_bb):
    for i2 in range(N_hh):
        Q[i1,i2], u_sol_0, u_sol_0_coarse, I_xx[i1,i2] = ChannelFlow(N_eta=N[0], N_xi=N[0], bb=bb[i1], hh=hh[i2], ll=ll, N_base=N[0])
Q = Q.T
I_xx = I_xx.T
bb_grid, hh_grid = np.meshgrid(bb, hh)
bb_grid = bb_grid
hh_grid = hh_grid

# contour_ticks = np.arange(0.00, 0.22, 0.02)
# color_ticks = np.arange(0.02, 0.22, 0.02)
fig, ax = plt.subplots()
# cs = ax.contour(x_grid, y_grid, u_sol, levels=color_ticks, colors='red', linestyles='dashed')
# plt.clabel(cs, inline=True, fontsize=10)
cp = ax.contourf(bb_grid, hh_grid, I_xx) # levels=contour_ticks)
cbar = fig.colorbar(cp)
plt.savefig("Fig1G_I.png", dpi=600)
ax.set_xlabel("b")
ax.set_ylabel("h")
ax.set_title("Finite difference solution for I")
plt.show()
