function [res] = jovanovic(parameters)


    % This function solves the Jovanovic (1979) model using value function iteration
    % Inputs:
    %   parameters: a structure with the following fields

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
    [m_grid, G_m]      = discretizeAR1_Tauchen(mu, 0, sigma_m,n_m, 3);          % m distribution
    m_grid = m_grid';
    [theta_0, F_theta] = discretizeAR1_Tauchen(0, 0, sigma_1,n_theta, 3);       % theta distribution (centered at 0)
    theta_0 = theta_0';

    % Keep only the firt row of the distribution
    G_m = G_m(1, :);
    F_theta = F_theta(1, :);

    % Initial guess for the value function
    V_m = zeros(1, n_m);
    theta_hat = zeros(1, n_m);
    m_hat = zeros(1, n_m);

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
            % Calculate the threshold theta_hat the lowest theta such that: (theta / (1 - beta) ≥ beta * Q)
            theta_pol = theta( find( theta / (1 - beta) > beta * Q, 1,  'first') );
            if numel(theta_pol) == 1
                theta_hat(i) = theta_pol;
            else
                theta_hat(i) = NaN;
            end
            % Calculate the threshold m_hat the lowest m such that: m_grid(i) + beta * max( int, beta * Q ) ≥ beta * Q
            m_pol = m_grid( find( m_grid(i) + beta * max( int, beta * Q ) >= beta * Q, 1,  'first') );
            if numel(m_pol) == 1
                m_hat(i) = m_pol;
            else
                m_hat(i) = NaN;
            end
        end

        % compute the distance between the old and the new value function
        dist = max(abs(V_m_new - V_m));

        % Update the value function
        V_m = V_m_new;
        
        % Update the iteration counter
        iter = iter + 1;
    end
    
    % Create a results structure
    res.V_m = V_m;
    res.theta_hat = theta_hat;
    res.m_hat = m_hat;

end  %jovanovic