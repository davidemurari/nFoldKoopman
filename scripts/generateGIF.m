function generateGIF(solution, lengths, n, dt, gif_filename, energy)
    % Generates a GIF showing the motion of a chain of connected pendula.
    %
    % Parameters:
    %   solution - (2n x N) matrix, where the first n rows are angles (theta)
    %              and the last n rows are momenta (p).
    %   lengths - Vector of pendulum lengths [l1, l2, ..., ln].
    %   n - Number of pendula.
    %   dt - Time step between frames.
    %   gif_filename - Name of the output GIF file.
    %   energy - Vector of Hamiltonian energy values (1 x N).

    % Extract angles (first n rows of the solution)
    theta = solution(1:n, :);
    N = size(theta, 2); % Number of time steps

    % Compute the positions of the pendula at each time step
    [x, y] = computePositions(theta, lengths);

    % Initialize figure
    figure;
    hold on;
    axis equal;
    xlim([-sum(lengths), sum(lengths)]); % Set x-axis limits
    ylim([-sum(lengths), sum(lengths)]); % Set y-axis limits
    xlabel('X');
    ylabel('Y');
    title('Chain Pendulum Simulation');

    % Create the GIF
    for i = 1:N
        % Clear the plot
        cla;

        % Plot the pendulum arms and masses
        x_current = [0; x(:, i)]; % Include the origin
        y_current = [0; y(:, i)]; % Include the origin

        % Plot the pendulum arms
        plot(x_current, y_current, 'r-', 'LineWidth', 2);

        % Plot the pendulum masses
        plot(x(1:n, i), y(1:n, i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

        % Update the title with the current time and energy
        current_time = (i - 1) * dt; % Calculate time from frame index
        current_energy = energy(i); % Get the energy at the current frame
        title(sprintf('Time: %.2f s | Energy: %.4f', current_time, current_energy));

        % Capture the frame
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        % Write to the GIF file
        if i == 1
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', dt);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', dt);
        end
    end

    disp(['GIF saved as ', gif_filename]);
end
