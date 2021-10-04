# Massachusetts Institute of Technology
# Department of Mechanical Engineering
author = "Zi Hao Foo"
date = "Oct 2021"

import numpy as np
import scipy as sp
import scipy.linalg as linalg
import matplotlib.pyplot as plt
import pandas as pd

# Supporting functions
def func_1(x_coords):
    """ returns f(x) = (6x + 2x^2) * exp(x) & u(x) = 2x * (1-x) * exp(x)"""
    f_hat = (6.0 * x_coords + 2.0 * (x_coords**2.0)) * np.exp(x_coords)
    u_actual = 2.0 * x_coords * (1.0 - x_coords) * np.exp(x_coords)
    return f_hat, u_actual

def func_2(x_coords):
    """ returns f(x) & u(x)"""
    f_hat = np.zeros(len(x_coords))
    u_actual = np.zeros(len(x_coords))
    for j1 in range(len(x_coords)):
        if x_coords[j1] < 0.5:
            f_hat[j1] = 2.0
            u_actual[j1] = 0.25 * x_coords[j1] * (3 - 4.0 * x_coords[j1]) 
        elif x_coords[j1] < 1.0:
            f_hat[j1] = 0.0
            u_actual[j1] = 0.25 * (1 - x_coords[j1]) 
    return f_hat, u_actual

def func_3(x_coords):
    """ returns f(x) = 3x & u(x) = -0.5x^3 + 0.5x"""
    f_hat = 3.0 * x_coords
    u_actual = -0.5 * (x_coords**3.0) + 0.5 * (x_coords)
    return f_hat, u_actual

# Main script
del_x_points = 8            # number of delta x 
del_x = np.array([0.1 * (0.5**i1) for i1 in range(del_x_points)])      # Array of delta_x, starting from 1/10 to 1/1280
num_points = (1.0 / del_x) - 1
num_points = num_points.astype(int)                         # Type cast to int
# print(num_points)
err_l1_norm = np.zeros(del_x_points)
err_l2_norm = np.zeros(del_x_points)
err_inf_norm = np.zeros(del_x_points)

for k1 in range(len(del_x)):
    A = 2.0*np.diag(np.ones(num_points[k1])) - np.diag(np.ones(num_points[k1]-1), k=1) - np.diag(np.ones(num_points[k1]-1), k=-1)       # Creates toeplitz matrix A
    x_coords = np.arange(0+del_x[k1], 1, del_x[k1])        # Creates x coordinates of the nodes
    # print("x coordinates:", x_coords)
    f_hat, u_actual = func_3(x_coords)                               # Discrete values of function f at x_coordinates
    # print("f hat values:", f_hat)
    # print("A shape:", A.shape)
    u_hat = linalg.solve(A, f_hat) * (del_x[k1]**2.0)
    # print("u hat valus:", u_hat)

    # Visualize FD approximation
    plt.plot(x_coords, u_hat, 'o')
    plt.plot(x_coords, u_actual, 'k-')
    plt.ylabel("u(x)")
    plt.xlabel("x")
    plt.show()
    
    err_vec = np.abs(u_hat - u_actual)
    err_l1_norm[k1] = del_x[k1] * linalg.norm(err_vec, 1)
    err_l2_norm[k1] = del_x[k1] * linalg.norm(err_vec, 2)
    err_inf_norm[k1] = del_x[k1] * linalg.norm(err_vec, np.inf)

alp_l1 = np.zeros(len(del_x)-1)
alp_l2 = np.zeros(len(del_x)-1)
alp_inf = np.zeros(len(del_x)-1)

for k2 in range(len(del_x)-1):
    alp_l1[k2] = np.log10(err_l1_norm[k2+1] / err_l1_norm[k2]) / np.log10(del_x[k2+1] / del_x[k2])
    alp_l2[k2] = np.log10(err_l2_norm[k2+1] / err_l2_norm[k2]) / np.log10(del_x[k2+1] / del_x[k2])
    alp_inf[k2] = np.log10(err_inf_norm[k2+1] / err_inf_norm[k2]) / np.log10(del_x[k2+1] / del_x[k2])

iter_coord = [j1 for j1 in range(len(del_x)-1)]
plt.plot(iter_coord, alp_l1, 'k', label = "alpha_1")
plt.plot(iter_coord, alp_l2, 'r', label = "alpha_2")
plt.plot(iter_coord, alp_inf, 'b', label = "alpha_inf")
plt.legend(loc="best")
plt.ylabel("Alpha")
plt.xlabel("Successive trials")
plt.show()