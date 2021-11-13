%% Zi Hao Foo
% 6.339 Project 2 (Finite Element Method)

clear; clc;

load grids

mu_0 = [0.4, 0.6, 0.8, 1.2, 0.1]; 
mesh = medium;
mu = mu_0;
        
        
%% Finite Element Solver
[u, Troot] = ThermalFin(mesh, mu);

%% Plotting
plotsolution(mesh, u); 