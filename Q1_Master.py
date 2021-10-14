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
from Q1_Sub import ChannelFlow

N = np.array([11, 21, 41, 81, 161, 321])
bb = 0.5        # in x domain
hh = 1.0        # in y domain
ll = 3.0
Q = np.zeros(len(N), dtype=float)
u_sol = np.zeros((len(N), max(N**2)))

for i1 in range(len(N)):
    Q[i1], u_sol = ChannelFlow(N_eta=N[i1], N_xi=N[i1], bb=bb, hh=hh, ll=ll)

y_N = N**(-2.0)
fig, ax = plt.subplots()
ax.loglog(N, Q, color='red')
# ax.loglog(N, y_N, color='black')
plt.show()

# print(u_sol)