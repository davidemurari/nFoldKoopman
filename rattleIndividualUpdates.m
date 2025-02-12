function [X, P, H] = rattle_simple_pendulum(dt, T, x0, p0, m, l, g)
    % dt: Time step
    % T: Total simulation time
    % x0: Initial positions (2N x 1)
    % p0: Initial momenta (2N x 1)
    % m: Mass vector (N x 1)
    % l: Length vector (N x 1)
    % g: Gravity constant
    
    N = floor(length(x0)/2); %number of pendula
    num_steps = floor(T / dt);
    e2 = [0;1];
    
    % Initialize storage
    X = zeros(2 * N, num_steps+1);
    P = zeros(2 * N, num_steps+1);
    H = zeros(1, num_steps+1);
    
    vec_e2 = kron(ones(N,1),e2);

    energy = @(u,v) sum(v.*v)/(2*m) + m*g*sum(vec_e2.*u);

    X(:,1) = x0;
    P(:,1) = p0;
    H(1) = energy(x0,p0);

    % Time-stepping with RATTLE
    for n = 1:num_steps
        x = X(:,n);
        p = P(:,n);

        % disp(['n: ' num2str(100*n/num_steps) '%'])

        %p_half_bar = @(lam) p - dt/2 * x * lam - dt/2 * g * m * e2;
        %x_new = @(lam) x + dt * p_half_bar(lam)/m;

        %x_update = x + dt * p/m - dt^2/2 * g * vec_e2;
        p_mid = p - dt/2 * m * g * vec_e2;
        x_update = x + dt/m * p_mid;
        tol = 1e-14;

        constraints_q = zeros(N,1);
        for i = 1:N
            if i==1
                constraints_q(1) = l^2-norm(x_update(1:2))^2;
            else
                constraints_q(i) = l^2-norm(x_update(2*i-1:2*i) - x_update(2*i-3:2*i-2))^2;
            end
        end

        it = 1;
        max_it = 1000;
        while norm(constraints_q,'inf')>tol && it<max_it
            for i = 1:N

                if i==1 && abs(constraints_q(1))>tol

                    lambda_1 = m/(dt^2) * constraints_q(1) / sum(x_update(1:2).*x(1:2));
                    x_update(1:2) = x_update(1:2) + dt^2/(2*m)*x(1:2)*lambda_1;
                    p_mid(1:2) = p_mid(1:2) + dt/2*x(1:2)*lambda_1;
                    constraints_q(1) = l^2-norm(x_update(1:2))^2;
                    
                end
            
                if i>1 && abs(constraints_q(i))>tol
                    lambda_i = m/(2*dt^2) * constraints_q(i) / sum((x_update(2*i-3:2*i-2)-x_update(2*i-1:2*i)).*(x(2*i-3:2*i-2)-x(2*i-1:2*i)));                    
                 
                    x_update(2*i-1:2*i) = x_update(2*i-1:2*i) + dt^2/(2*m)*(x(2*i-1:2*i)-x(2*i-3:2*i-2))*lambda_i;
                    x_update(2*i-3:2*i-2) = x_update(2*i-3:2*i-2) - dt^2/(2*m)*(x(2*i-1:2*i)-x(2*i-3:2*i-2))*lambda_i;
                    
                    p_mid(2*i-1:2*i) = p_mid(2*i-1:2*i) + dt/2*(x(2*i-1:2*i)-x(2*i-3:2*i-2))*lambda_i;
                    p_mid(2*i-3:2*i-2) = p_mid(2*i-3:2*i-2) - dt/2*(x(2*i-1:2*i)-x(2*i-3:2*i-2))*lambda_i;

                    constraints_q(i) = l^2-norm(x_update(2*i-3:2*i-2)-x_update(2*i-1:2*i))^2;
                end
            end

            it = it + 1;
        end

        %p_mid = m*(x_update - x) / dt; %p - dt/2 * m * g * vec_e2 - dt/2 * G' * lambda;
        p_update = p_mid - dt/2 * m * g * vec_e2;

        constraints_p = zeros(N,1);
        for i = 1:N
            if i==1
                constraints_p(1) = sum(p_update(1:2).*x_update(1:2));
            else
                constraints_p(i) = sum((p_update(2*i-3:2*i-2)-p_update(2*i-1:2*i)).*(x_update(2*i-3:2*i-2)-x_update(2*i-1:2*i)));
            end
        end

        it = 1;
        max_it = 1000;
        while norm(constraints_p,'inf')>tol && it<max_it

            for i = 1:N
                
                if i==1 && abs(constraints_p(i))>tol
                    p_mid1 = p_update(1:2);
                    q_new_1 = x_update(1:2);
                    mu_1 = m/dt * constraints_p(1) / l^2;
                    p_mid1 = p_mid1 - dt/2 * mu_1 * q_new_1;
                    constraints_p(i) = sum(p_mid1.*q_new_1);
                    p_update(1:2) = p_mid1;
                end

                if i>1 && abs(constraints_p(i))>tol
                    pi = p_update(2*i-1:2*i);
                    pi1 = p_update(2*i-3:2*i-2);
                    qi = x_update(2*i-1:2*i); %old points
                    qi1 = x_update(2*i-3:2*i-2); %old points

                    mu_i = m/dt * constraints_p(i) / l^2;
                    pi = pi - dt/2 * (qi-qi1) * mu_i;
                    pi1 = pi1 + dt/2 * (qi-qi1) * mu_i;

                    constraints_p(i) = sum((pi1-pi).*(qi1-qi));
                    p_update(2*i-1:2*i) = pi;
                    p_update(2*i-3:2*i-2) = pi1;
                end
            end
            it = it + 1;
        end

        % disp("==== CHECK CONSTRAINTS =====")
        % disp(['ref length is ' num2str(l)])
        % for i = 1:N
        %     if i==1
        %         disp(['i: ' num2str(i)])
        %         disp(['length ' num2str(norm(x_update(1:2))^2-l^2)]);
        %         disp(['tangent ' num2str(sum(x_update(1:2).*p_update(1:2)))]);
        %     else
        %         disp(['i: ' num2str(i)])
        %         disp(['computed length ' num2str(norm(x_update(2*i-3:2*i-2) - x_update(2*i-1:2*i)))])
        %         disp(['length ' num2str(norm(x_update(2*i-3:2*i-2) - x_update(2*i-1:2*i))^2-l^2)]);
        %         disp(['tangent ' num2str(sum((p_update(2*i-3:2*i-2)-p_update(2*i-1:2*i)).*(x_update(2*i-3:2*i-2)-x_update(2*i-1:2*i))))]);
        %     end
        % end
        % disp("==== END OF THE CHECK ====")

        % Store results
        X(:,n+1) = x_update;
        P(:,n+1) = p_update;
        H(n+1) = energy(x_update,p_update);
    end
end

% Simulation setup
N = 50;
L = 1;
M = 1;
g = 9.81;
dt = 5e-4;
final = 1.;
T = 0.2; %floor(final/dt)*dt;
l = L / N; %* ones(N, 1);
m = M / N; %* ones(N, 1);

number_samples = 100;

X_data = zeros(2*N,number_samples);
Y_data = zeros(2*N,number_samples);

for it = 1:number_samples

    qq = 2*rand(N,2)-1;
    qq = l * (qq ./ vecnorm(qq,2,2)); %l2 norm row-wise
    qq = cumsum(qq,1);
    qq = reshape(qq.',1,[])'; %write as a single column stacking the rows
    
    pp = zeros(2*N,1);
    J = randn(2,2);
    J = (J-J')/2;
    pp(1:2) = J * qq(1:2);
    for i = 1:N-1
        J = randn(2,2);
        J = (J-J')/2;
        v = (qq(2*i+1:2*i+2) - qq(2*i-1:2*i));
    
        pp(2*i+1:2*i+2) = J * v + pp(2*i-1:2*i);
    end
    % Initial conditions (hanging vertically at rest)
    %x0 = [0; -l(1); l(2)*sin(pi/4); -l(1)-l(2)*cos(pi/4)];
    %x0 = kron(cumsum(l*ones(N,1)),[cos(pi/12);sin(pi/12)]);
    %x0 = l * [sin(pi/6);cos(pi/6); 0 + sin(pi/6); 1 + cos(pi/6)];
    
    %p0 = zeros(2*N, 1);
    
    % e2 = [0;1];
    % vec_e2 = kron(ones(N,1),e2);
    % energy = @(u,v) sum(v.*v)/(2*m) + m*g*sum(vec_e2.*u);
    % E0 = energy(x0,p0);
    % v_star = sqrt(2*E0/m);
    % dt = l / (5*v_star);
    % disp(['dt' num2str(dt)])
    
    % disp("==== CHECK CONSTRAINTS =====")
    % disp(['ref length is ' num2str(l)])
    % for i = 1:N
    %     if i==1
    %         disp(['i: ' num2str(i)])
    %         disp(['length ' num2str(norm(qq(1:2))^2-l^2)]);
    %         disp(['tangent ' num2str(sum(qq(1:2).*pp(1:2)))]);
    %     else
    %         disp(['i: ' num2str(i)])
    %         disp(['computed length ' num2str(norm(qq(2*i-3:2*i-2) - qq(2*i-1:2*i)))])
    %         disp(['length ' num2str(norm(qq(2*i-3:2*i-2) - qq(2*i-1:2*i))^2-l^2)]);
    %         disp(['tangent ' num2str(sum((pp(2*i-3:2*i-2)-pp(2*i-1:2*i)).*(qq(2*i-3:2*i-2)-qq(2*i-1:2*i))))]);
    %     end
    % end
    
    
    % Run simulation
    tic
    [X, P, H] = rattle_simple_pendulum(dt, T, qq, pp, m, l, g);
    toc
    disp(['progress bar: ' num2str(100*it/number_samples) '%'])
    X_data(:,it) = X(:,1);
    Y_data(:,it) = X(:,end);
end

filename = sprintf('X_data_%d_pendula_%d_samples.csv', N,number_samples);
writematrix(X_data, filename);
filename = sprintf('Y_data_%d_pendula_%d_samples.csv', N,number_samples);
writematrix(Y_data, filename);

% % Plot results
% figure;
% plot(squeeze(X(1,:)), squeeze(X(2,:)), 'ro', 'LineWidth', 2); hold on;
% xlabel('x'); ylabel('y'); title('Trajectory');
% grid on;
% axis equal
% hold off;
% 
% time_array = 0:dt:T;
% 
% figure;
% plot(time_array',H,'k','LineWidth', 2);
% xlabel('t'); ylabel('H'); title('Hamiltonian energy');
% grid on;
% hold off;
% 
% figure;
% semilogy(time_array',abs(H-H(1)),'k','LineWidth', 2);
% xlabel('t'); ylabel('H'); title('Hamiltonian energy variation');
% grid on;
% hold off;
% 
% toc
% 
% % Create figure
% fig = figure;
% filename = 'double_pendula.gif';
% 
% % Ensure the figure is visible (optional)
% set(fig, 'Visible', 'on');
% 
% for n = 1:10:length(time_array)
%     % Clear the figure
%     clf;
% 
%     % Extract the positions of each pendulum node
%     % (X is 2N x (num_steps+1), so each rod end is (X(2i-1,n), X(2i,n)))
%     hold on;
%     for i = 1:N
%         x_coord = X(2*i - 1, n);
%         y_coord = X(2*i,     n);
% 
%         % Plot each pendulum node
%         plot(x_coord, y_coord, 'ro', 'MarkerSize', 2, 'MarkerFaceColor', 'r');
% 
%         % You might also want to connect them with lines (i.e. rods):
%         %   If the first rod is from (0,0) to the first node:
%         if i == 1
%             plot([0, x_coord], [0, y_coord], 'b-','LineWidth',1.5);
%         else
%             % The i-th rod is from node i-1 to node i
%             x_im1 = X(2*(i-1)-1, n);
%             y_im1 = X(2*(i-1),   n);
%             plot([x_im1, x_coord], [y_im1, y_coord], 'b-','LineWidth',1.5);
%         end
%     end
%     hold off;
% 
%     % Set axis limits & aspect ratio as desired
%     axis equal;
%     axis([-2 2 -2 2]);  % or something suitable to your system
% 
%     % Title with current energy
%     title(sprintf('Time %.2f,  Energy = %.4f', time_array(n), H(n)));
% 
%     % Force drawing & capture the frame
%     drawnow;
%     frame = getframe(fig);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,256);
% 
%     % Write/append to GIF
%     if n == 1
%         % First frame initializes the GIF file
%         imwrite(A,map,filename,'gif','LoopCount',inf,'DelayTime',0.05);
%     else
%         % Subsequent frames are appended
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
%     end
% end
