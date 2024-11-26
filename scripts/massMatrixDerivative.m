function dM_dtheta_j = massMatrixDerivative(j, n, lengths, masses, theta)
    % Computes the derivative of the mass matrix M(theta) with respect to theta_j.
    %
    % Parameters:
    %   j - Index of the angle with respect to which the derivative is taken (1-based index).
    %   n - Number of pendulums.
    %   lengths - Vector of lengths of the pendulums [l1, l2, ..., ln].
    %   masses - Vector of masses of the pendulums [m1, m2, ..., mn].
    %   theta - Vector of angles of the pendulums [theta1, theta2, ..., thetan].
    %
    % Returns:
    %   dM_dtheta_j - The derivative dM/dtheta_j as an (n, n) matrix.

    % Initialize the derivative matrix with zeros
    dM_dtheta_j = zeros(n, n);

    % Loop through rows and columns
    for i = 1:n
        for k = 1:n
            % Compute the sum of masses from max(i, k) to n
            mass_sum = sum(masses(max(i, k):n));

            % Compute the derivative terms
            if i == j
                dM_dtheta_j(i, k) = -lengths(i) * lengths(k) * sin(theta(i) - theta(k)) * mass_sum;
            end
            if j == k
                dM_dtheta_j(i, k) = lengths(i) * lengths(k) * sin(theta(i) - theta(k)) * mass_sum;
            end
        end
    end
end