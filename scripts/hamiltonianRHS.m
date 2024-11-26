function dydt = hamiltonianRHS(state, n, lengths, masses, g)
    % Computes the right-hand side of the Hamiltonian equations.
    %
    % Parameters:
    %   t - Time (not explicitly used, included for compatibility with ODE solvers).
    %   state - Current state vector [theta1, ..., thetan, p1, ..., pn].
    %   n - Number of pendulums.
    %   lengths - Vector of lengths of the pendulums.
    %   masses - Vector of masses of the pendulums.
    %   g - Gravitational acceleration.
    %
    % Returns:
    %   dydt - Derivatives [dtheta1, ..., dthetan, dp1, ..., dpn].

    % Split state vector into angles (theta) and momenta (p)
    theta = state(1:n);
    p = state(n+1:end);

    % Assemble the mass matrix M
    M = massMatrix(n, lengths, masses, theta);
    
    % Solve M * theta_dot = p for theta_dot
    theta_dot = M \ p;

    % Compute the potential energy gradient dV/dtheta
    dV_dtheta = zeros(n, 1);
    for i = 1:n
        % Sum of masses from index i to n
        mgj = sum(masses(i:n)) * g;
        dV_dtheta(i) = mgj * lengths(i) * sin(theta(i));
    end

    % Compute the derivative of kinetic energy dT/dtheta
    dT_dtheta = zeros(n, 1);
    for i = 1:n
        % Compute the derivative of the mass matrix with respect to theta_i
        dM_dtheta_i = massMatrixDerivative(i, n, lengths, masses, theta);

        % Compute y = M^-T * (dM/dtheta_i * M^-1 * p)
        y = M \ (dM_dtheta_i * theta_dot); % Equivalent to solve(M^T, ...)

        % Compute the contribution to dT_dtheta
        dT_dtheta(i) = -0.5 * (p' * y);
    end

    % Compute p_dot = -dV/dtheta - dT/dtheta
    p_dot = -dV_dtheta - dT_dtheta;

    % Concatenate derivatives for theta and p
    dydt = [theta_dot; p_dot];
end
