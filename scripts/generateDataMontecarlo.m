function [X,Y] = generateDataMontecarlo(n_samples, n_pendula, n_steps, h, plb, pub)
    
    
    % System parameters
    L = 1.0;  % Total length of the chain of pendulums
    m = 1.0;  % Total mass of the chain of pendulums
    lengths = L / n_pendula * ones(n_pendula, 1); % Uniform lengths
    masses = m / n_pendula * ones(n_pendula, 1); % Uniform masses
    g = 9.81;  % Gravitational acceleration
    
    % Time span
    t0 = 0; % Initial time
    tf = n_steps * h; % Final time
    
    % Define the Hamiltonian RHS and function
    XH = @(z) hamiltonianRHS(z, n_pendula, lengths, masses, g);
    H = @(z) hamiltonianFunction(n_pendula, lengths, masses, g, z(1:n_pendula), z(n_pendula+1:end));
    
    data = zeros(n_samples, n_steps + 1, 2 * n_pendula);
    X = zeros(2*n_pendula, n_samples*n_steps);
    Y = zeros(2*n_pendula, n_samples*n_steps);

    hWaitBar = waitbar(0, 'Processing steps...', 'Name', 'Progress', ...
                   'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');

    setappdata(hWaitBar, 'canceling', 0); % Add a canceling flag

    
    for i = 1:n_samples
        is_not_converged = 1;
        while is_not_converged == 1
            ic = sampleICs(1, n_pendula, plb, pub); % Sample n_samples ICs
            [cc, ~, is_not_converged] = implicit_midpoint_solver(XH, H, ic', t0, tf, n_steps);
            if is_not_converged==0
                data(i, :, :) = cc';
            end
        end
        waitbar(i / n_samples, hWaitBar, sprintf('Generated trajectory %d of %d', i, n_samples));
    end

    delete(hWaitBar);

    for i = 1:n_samples
        X(:,(i-1)*n_steps + 1:i*n_steps) = squeeze(data(i,1:n_steps,:))';
        Y(:,(i-1)*n_steps + 1:i*n_steps) = squeeze(data(i,2:n_steps+1,:))';
    end
end

    
