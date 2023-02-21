function [res] = jovanovic(parameters)


    % This function solves the Jovanovic (1979) model using value function iteration
    % Inputs:
    %   parameters: a structure with the following fields
    %       beta: discount factor
    % Output:
    %   res: a structure with the following fields
    %       V_m: value function
    %       theta_hat: threshold theta
    %       m_hat: threshold m
    

    % Unpack the parameters
    beta = parameters.beta;
    mu = parameters.mu;
    sigma_0 = parameters.sigma_0;
    sigma_u = parameters.sigma_u;
    n_m = parameters.n_m;
    n_theta = parameters.n_theta;
    tol = parameters.tol;
    n_iter = parameters.n_iter;
    
    
    % Distributions
    sigma_m = ( sigma_0^4 ) / ( sigma_0^2 + sigma_u^2 );                        % m distribution
    sigma_1 = ( sigma_0^2 ) / ( sigma_0^2 + sigma_u^2 ) * sigma_u^2;            % theta distribution

    % Grids
    m_grid = linspace(mu - 3 * sigma_m, mu + 3 * sigma_m, n_m);
    theta_0 = linspace(- 3 * sigma_1, 3 * sigma_1, n_theta);
    G_m = normpdf(m_grid, mu, sigma_m); G_m = G_m / sum(G_m);
    F_theta = normpdf(theta_0, 0, sigma_1); F_theta = F_theta / sum(F_theta);

    % Add theta

    % Value function initial guess
    V_m = zeros(1, n_m);
    
    %% Solving the model
    % Initialize the error and the iteration counter
    dist = inf;
    iter = 0;

    % Iterate until the error is smaller than the tolerance
    while (dist > tol) && (iter < n_iter)

        % Calculate Q
        Q = V_m * G_m';

        % Calculate the value function

        % Allocate memory for the value function
        V_m_new = zeros(1, n_m);
        % Iteate over the m grid
        for i = 1:n_m
            % theta distribution just add m to the values
            theta = theta_0 + m_grid(i);
            % Calculat the integral over theta when theta > Q
            int_over = ( (theta / (1 - beta)) .* (theta / (1 - beta) > beta * Q) ) * F_theta';
            % Calculate the integral over theta when theta < Q 
            int_under = ( beta * Q ) * ( (theta / (1 - beta) < beta * Q) * F_theta' );
            % Add the two integrals
            int = int_over + int_under;
            % Calculate the value function
            V_m_new(i) = max(m_grid(i) + beta * max( int, beta * Q ), beta * Q);
        end

        % compute the distance between the old and the new value function
        dist = max(abs(V_m_new - V_m));

        % Update the value function
        V_m = V_m_new;
        
        % Update the iteration counter
        iter = iter + 1;
    end

    % Calculate the threshold theta_hat the lowest theta such that: (theta / (1 - beta) ≥ beta * Q)
    theta_hat = theta( find( theta / (1 - beta) > beta * Q, 1,  'first') );

    % Calculate the threshold m_hat the lowest m such that: m_grid(i) + beta * max( int, beta * Q ) ≥ beta * Q
    m_hat = m_grid( find( m_grid(i) + beta * max( int, beta * Q ) >= beta * Q, 1,  'first') );
    
    % Create a results structure
    res.V_m = V_m;
    res.theta_hat = theta_hat;
    res.m_hat = m_hat;

end  %jovanovic