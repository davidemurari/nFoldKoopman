function [x,y] = computePositions(theta,lengths)
    [n, num_steps] = size(theta);
    x = zeros(n, num_steps);
    y = zeros(n, num_steps);
    for i = 1:n
        x(i, :) = sum(lengths(1:i) .* sin(theta(1:i, :)), 1);
        y(i, :) = -sum(lengths(1:i) .* cos(theta(1:i, :)), 1);
    end
end