function [z_sol, t_sol, energy] = ode45_solver(XH, H, z0, t0, tf, N)
    % Inputs:
    % XH - Hamiltonian vector field
    % H - function handle for the Hamiltonian H(q, p)
    % z0 - initial condition, a 2*d vector [q0; p0]
    % t0, tf - start and end times
    % N - number of output points for uniform sampling
    
    % Ensure z0 has even length
    assert(mod(length(z0), 2) == 0, 'z0 must have an even number of elements');
    d = length(z0) / 2; % Dimension of q or p

    % Define the ODE system
    ode_system = @(t, z) XH(z);

    % Solve the system using ode45
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
    [t_sol, z_sol] = ode45(ode_system, [t0, tf], z0, options);

    z_sol = z_sol';

    % Interpolate solution to uniform time vector
    %z = interp1(t_sol, z_sol, t_out)';

    % Compute energy at each time step
    energy = zeros(1, length(t_sol));
    for n = 1:length(t_sol)
        energy(n) = H(z_sol(:,n));
    end
end
