function [z,energy,is_not_converged]=implicit_midpoint_solver(XH, H, z0, t0, tf, N)
    % Inputs:
    % XH - Hamiltonian vector field
    % H - function handle for the Hamiltonian H(q, p)
    % z0 - initial condition, a 2*d vector [q0; p0]
    % t0, tf - start and end times
    % N - number of time steps
    
    % Extract dimensions
    assert(mod(length(z0), 2) == 0, 'z0 must have an even number of elements');
    d = length(z0) / 2;      % Dimension of q or p
    
    is_not_converged = 0;

    % Time step and time vector
    h = (tf - t0) / (N);
    
    % Initialize arrays to store solutions
    z = zeros(2*d, N+1);
    energy = zeros(1, N+1);
    z(:, 1) = z0;
    
    % Initial energy
    energy(1) = H(z0);
    
    % Options for fsolve
    options = optimoptions('fsolve', 'TolFun', 1e-8, 'TolX', 1e-8, 'Display', 'none');
    substeps = floor(h/0.01); %Number of substeps done to generate a data point
    % Implicit midpoint method loop
    for n = 1:N
        % Define the nonlinear system to solve
        supp = z(:,n);
        for it = 1:substeps
            system = @(mid) mid - supp - h/substeps * XH(0.5*(mid+supp));
            [supp, ~, exitflag] = fsolve(system, supp, options);
            if exitflag<=0
                is_not_converged = 1;
                break
            end
        end
        z(:,n+1) = supp;
        energy(n+1) = H(z(:,n+1));
    end
   
end
