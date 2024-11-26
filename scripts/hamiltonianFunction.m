function H=hamiltonianFunction(n, lengths, masses, g, theta, p)
    [x, y] = computePositions(theta, lengths);

    % Assemble the mass matrix
    M = massMatrix(n, lengths, masses, theta);

    % Compute potential energy (V)
    V = 0;
    for i = 1:n
        V = V + masses(i) * g * y(i);
    end

    % Compute kinetic energy (T)
    T = 0.5 * (p' * (M \ p)); % Equivalent to np.dot(p, solve(M, p)) in Python

    % Hamiltonian is the sum of kinetic and potential energy
    H = T + V;
end