function [params] = createParams(n_m, n_theta, sigma_0, sigma_u, mu, beta, n_iter, tol);


    % Create a structure of parameters for the simulation
    
    params.n_m = n_m;                                   % number of m-grid points
    params.n_theta = n_theta;                           % number of theta-grid points
    params.sigma_0 = sigma_0;                           % standard deviation of -----
    params.sigma_u = sigma_u;                           % standard deviation of -----
    params.mu = mu;                                     % mean of the prior
    params.beta = beta;                                 % precision of the noise
    params.n_iter = n_iter;                             % number of iterations
    params.tol = tol;                                   % tolerance for the stopping criterion


end  %createParams