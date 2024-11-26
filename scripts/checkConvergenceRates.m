% Example: Simple pendulum
XH = @(t, z) [z(2); -sin(z(1))]; % Hamiltonian system (vector field)
H = @(z) 0.5 * z(2)^2 - cos(z(1)); % Energy function

% Initial conditions
z0 = [5/4*pi; 0.5]; % Initial state
t0 = 0;      % Initial time
tf = 1;     % Final time
Ns = 2.^[8,9,10]; % Number of steps to test

% Generate reference solution using ode45
ref_sol = ode45(XH, [t0 tf], z0);
ref_times = linspace(t0, tf, 1000);
ref_values = deval(ref_sol, ref_times);

solvers = {
    @(XH, H, z0, t0, tf, N) euler_heun_solver(XH, H, z0, t0, tf, N), 'Euler Heun';
    @(XH, H, z0, t0, tf, N) implicit_euler_solver(XH, H, z0, t0, tf, N), 'Implicit Euler';
    @(XH, H, z0, t0, tf, N) implicit_midpoint_solver(XH, H, z0, t0, tf, N), 'Implicit Midpoint';
    @(XH, H, z0, t0, tf, N) trapezoidal_solver(XH, H, z0, t0, tf, N), 'Trapezoidal Method';
    @(XH, H, z0, t0, tf, N) explicit_euler_solver(XH, H, z0, t0, tf, N), 'Explicit Euler';
    @(XH, H, z0, t0, tf, N) symplectic_euler_solver(XH, H, z0, t0, tf, N), 'Symplectic Euler'
};

for i = 1:size(solvers, 1)
    solver = solvers{i, 1};
    solver_name = solvers{i, 2};
    
    fprintf('Validating solver: %s\n', solver_name);

    % Compute errors for different step sizes
    errors = zeros(size(Ns));
    for j = 1:length(Ns)
        N = Ns(j);
        [z, energy] = solver(@(z) XH(0,z), H, z0, t0, tf, N); % Run the solver
        z = z';
        times = linspace(t0, tf, N + 1);            % Time grid used in solver
        ref_interp = interp1(ref_times', ref_values', times', 'linear'); % Interpolate ref solution
        
        errors(j) = vecnorm(ref_interp(end) - z(end), 2, 1);
    end

    figure;
    loglog(Ns, errors, '-o', 'LineWidth', 2, 'DisplayName', 'Numerical Error');
    hold on;

    % Compute the estimated orders of convergence
    orders = log(errors(1:end-1) ./ errors(2:end)) ./ log(Ns(2:end) ./ Ns(1:end-1));
    
    % Compute the average estimated order of convergence
    average_order = mean(orders);

    % Add reference lines with slopes 1, 2, 3, 4
    last_error = errors(end); % Error at the finest resolution
    last_N = Ns(end);         % Number of steps at the finest resolution
    for slope = 1:4
        ref_line = last_error * (Ns / last_N).^(-slope);
        loglog(Ns, ref_line, '--', 'LineWidth', 1.5, 'DisplayName', sprintf('Slope %d', slope));
    end
    
    xlabel('Number of Steps (N)');
    ylabel('Max Error');
    title(sprintf('Convergence of %s (Estimated Order: %.2f)', solver_name, average_order));
    legend('show', 'Location', 'best');
    grid on;
    hold off;
end
