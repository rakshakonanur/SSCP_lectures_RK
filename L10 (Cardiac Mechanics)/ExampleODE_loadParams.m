function [params,y0] = ExampleODE_loadParams() 
% ExampleODE_loadParams.m 
% Automatically generated by Netflux on 09-May-2025
 
% species parameters 
speciesNames = {'A','B','C','D','E',}; 
tau = [1, 1, 1, 1, 1, ]; 
ymax = [1, 1, 1, 1, 1, ];
 
% reaction parameters 
w = [1, 1, 1, 1, 1, 1, ]; 
n = [1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, ]; 
EC50 = [5.000000e-01, 5.000000e-01, 5.000000e-01, 5.000000e-01, 5.000000e-01, 5.000000e-01, ]; 
rpar = [w;n;EC50];
 
params = {rpar,tau,ymax,speciesNames};
 
y0 = [0, 0, 0, 0, 0, ]; 
