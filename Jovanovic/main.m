%{
==========================================================================================
Title: : Jovanovic Model of Job Matching
Author: Mitchell Valdes-Bobes @mitchv34
        Yeonggyu Yun
Date: 2023-02-16
Description: 
==========================================================================================
%}

% Parameters 
% Parameters
n_m = 100; % Number of m values
n_theta = 100; % Number of theta values
sigma_0 = 0.1;
sigma_u = 0.1;
mu = 0;
beta = 0.995;
n_iter = 1000;
tol = 1e-6;

% Create the parameter structure
params = createParams(n_m, n_theta, sigma_0, sigma_u, mu, beta, n_iter, tol);

% Solve the model
res = jovanovic(params);