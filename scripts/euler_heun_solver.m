function [z,energy]=euler_heun_solver(XH, H, z0, t0, tf, N)
    % Inputs:
    % XH - Hamiltonian vector field
    % H - function handle for the Hamiltonian H(q, p)
    % z0 - initial condition, a 2*d vector [q0; p0]
    % t0, tf - start and end times
    % N - number of time steps

    % Second order explicit method
    
    % Extract dimensions
    assert(mod(length(z0), 2) == 0, 'z0 must have an even number of elements');
    d = length(z0) / 2;      % Dimension of q or p
    
    % Time step and time vector
    h = (tf - t0) / (N);
    
    % Initialize arrays to store solutions
    z = zeros(2*d, N+1);
    energy = zeros(1, N+1);
    z(:, 1) = z0;
    
    % Initial energy
    energy(1) = H(z0);
    
    % Euler Heun method loop
    for n = 1:N
        zhat = z(:,n) + h * XH(z(:,n));
        z(:, n+1) = z(:,n) + h/2 * (XH(zhat) + XH(z(:,n)));
        energy(n+1) = H(z(:,n+1));
    end
   
end
