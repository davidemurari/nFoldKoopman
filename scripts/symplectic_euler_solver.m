function [z, energy] = symplectic_euler_solver(XH, H, z0, t0, tf, N)
    % Inputs:
    % XH - Hamiltonian vector field
    % H - function handle for the Hamiltonian H(q, p)
    % z0 - initial condition, a 2*d vector [q0; p0]
    % t0, tf - start and end times
    % N - number of time steps
    
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

    get_q = @(x) x(1:d);
    get_p = @(x) x(d+1:end);

    % Symplectic Euler method loop
    for n = 1:N
        q_n = z(1:d, n);        % Position vector q
        p_n = z(d+1:end, n);    % Momentum vector p
        
        % Update p implicitly
        p_eqn = @(p_next) p_next - (p_n + h * get_p(XH([q_n; p_next]))); 
        p_next = fsolve(p_eqn, p_n, optimoptions('fsolve', 'TolFun', 1e-12, 'TolX', 1e-12));
        
        % Update q explicitly
        q_next = q_n + h * get_q(XH([q_n; p_next]));
        
        % Store updated values
        z(:, n+1) = [q_next; p_next];
        energy(n+1) = H(z(:, n+1));
    end
end
