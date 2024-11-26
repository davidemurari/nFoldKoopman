function M = massMatrix(n, lengths, masses, theta)
    % Initialize an n x n zero matrix
    M = zeros(n, n);
    
    % Loop through rows and columns
    for i = 1:n
        for k = 1:n
            % Compute the sum of masses from max(i, k) to n
            mass_sum = sum(masses(max(i, k):n));
            
            % Compute the entry of the matrix
            M(i, k) = lengths(i) * lengths(k) * cos(theta(i) - theta(k)) * mass_sum;
        end
    end
end