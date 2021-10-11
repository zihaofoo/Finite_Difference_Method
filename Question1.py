# Massachusetts Institute of Technology
# Department of Mechanical Engineering
author = "Zi Hao Foo"
date = "Oct 2021"

import numpy as np
import scipy as sp
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import pandas as pd

# def ChannelFlow(N_xi, N_eta, bb, hh):

N_xi = 21
N_eta = 21 
bb = 0.5        # in x domain
hh = 1.0        # in y domain
ll = 3.0
a = np.sqrt(0.25 * ((ll - bb)**2) - (hh**2))
print(a)
# Inputs:
# N_xi : Number Of nodes in the xi direction.
# N_eta: Number of nodes int he eta direction
# bb : The length of the base
# hh : the height of the channel

# Outputs:
# Solution : A xi-eta matrix containing the solution
# Flowrate : The flowrate throught eh channel
# I_xx     : The moment of inertia of the channel

d_xi  = 1.0/(N_xi-1);            # delta xi                     
d_eta = 1.0/(N_eta-1);           # delta eta                      
NumNodes = N_xi * N_eta;   

# Constructing the Jacobian
J = np.zeros((N_xi, N_eta), dtype=float)
for i1 in range(N_xi): 
    for i2 in range(N_eta):  
        J[i1,i2] = hh * (0.5 * bb + a * i2)

# Initializing sparse matrices 'A' and 'Node' (coefficient value at i,j)
A = np.zeros((N_xi, N_eta), dtype=float)  
Node = np.zeros((N_xi, N_eta), dtype=float)                   

# Interior points







# return Solution, FlowRate, I_xx
