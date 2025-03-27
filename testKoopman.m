addpath('scripts');

% n_samples = 1000;
% n_pendula = 1;
% n_steps = 1;
h = 0.05; %to generate a data point the method will do floor(h/0.01) substeps
plb = -2; %Lower bound for the interval where momenta are sampled
pub = 2; %Upper bound for the interval where momenta are sampled
% 
% [X,Y] = generateDataMontecarlo(n_samples, n_pendula, n_steps, h, plb, pub);
% %X = sol_trpz(:,1:3000);
% %Y = sol_trpz(:,2:3001);

n_pendula = 1;
n_samples = 40000;
n_steps = 200;

%data = generateDataMontecarlo(n_samples, n_pendula, n_steps, h, plb, pub);
% return

filename = sprintf('X_data_%d_pendula_%d_samples_%d_steps.csv', n_pendula,n_samples,n_steps);
X_2d = readmatrix(filename);
filename = sprintf('Y_data_%d_pendula_%d_samples_%d_steps.csv', n_pendula,n_samples,n_steps);
Y_2d = readmatrix(filename);

X = reshape(X_2d, 4*n_pendula, n_samples, n_steps);
Y = reshape(Y_2d, 4*n_pendula, n_samples, n_steps);

% idx = max(abs(Y), [], 1) < 10 & ~any(isnan(Y), 1); %I am using this as a
% criterion for a converged solution % this will be hard with delay
% embedding so let's try and od without
% 
% X = X(:, idx);
% Y = Y(:, idx); %remove the non converged solutions

% n_samples = size(Y,2); 
% disp(n_samples)

%% MATRIX-BASED APPROACH

% %% Build matrix using kernel method
% [~,K,~,PX,PY] = kernel_ResDMD(X,Y,'type',"Gaussian",'N',500); % N = number of basis functions
% [~,LAM,W2] = eig(K,'vector');
% R = (sqrt(real(diag(W2'*L*W2)./diag(W2'*W2)-abs(LAM).^2))); % dual residual - warning, not L^2 residual

%% Build matrix using delay embedding
N = n_steps; % number of basis functions
M = n_samples; % number of data points
dim_sol = 4*n_pendula; %2*n_pendula;
PX = zeros(M,N*dim_sol);
PY = zeros(M,N*dim_sol);
% h = 1; % number of time steps for delay
% data = zeros(n_samples, n_steps + 1, 2 * n_pendula);
for jj= 1:N
    PX(:,(jj-1)*dim_sol+1:jj*dim_sol) = X(:,:,jj)';
    %PX(:,(jj-1)*dim_sol+1:jj*dim_sol) = squeeze(data(:,jj,:));
    PY(:,(jj-1)*dim_sol+1:jj*dim_sol) = Y(:,:,jj)';
    %PY(:,(jj-1)*dim_sol+1:jj*dim_sol) = squeeze(data(:,jj+1,:));
end

disp(["PX: " num2str(size(PX))])

% 
% % PX = X(:,1:end-1)';
% % PY = Y(:,1:end-1)';
% 
% disp(["Size PX: " size(PX)])
% disp(["Size PY: " size(PY)])
% 
% PX = X(:,1:end-1)';
% PY = Y(:,1:end-1)';
K = PX\PY;

disp(["Norm DMD approximation : " num2str(norm(K,2))])
[V,LAM] = eig(K,'vector');
R = vecnorm(PY*V-PX*V*diag(LAM))./vecnorm(PX*V); % L^2 residual


%% Plot eigenvalues and residuals

figure
scatter(real(LAM),imag(LAM),20,R,'filled','MarkerEdgeColor','k','LineWidth',0.01);
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-m')
axis equal
axis([-1.15,1.15,-1.15,1.15])
load('cmap.mat')
colormap(cmap2); colorbar
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title('Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

% do the same for mpEDMD
[~,mpV,mpD] = mpEDMDqr(PX,PY,1/M);
R2 = vecnorm(PY*mpV-PX*mpV*mpD)./vecnorm(PX*mpV);

figure
scatter(real(diag(mpD)),imag(diag(mpD)),20,R2,'filled','MarkerEdgeColor','k','LineWidth',0.01);
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'-m')
axis equal
axis([-1.15,1.15,-1.15,1.15])
load('cmap.mat')
colormap(cmap2); colorbar
xlabel('$\mathrm{Re}(\lambda)$','interpreter','latex','fontsize',18)
ylabel('$\mathrm{Im}(\lambda)$','interpreter','latex','fontsize',18)
title('Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

% figure
% loglog([0.001,1],[0.001,1],':k','linewidth',1)
% hold on
% loglog(sqrt(abs(abs(LAM).^2-1)),R,'k.','markersize',20)
% xlabel('$\sqrt{1-|\lambda|^2}$','interpreter','latex','fontsize',18)
% title('Dual Eigenpair Residuals','interpreter','latex','fontsize',18)
% ax=gca; ax.FontSize=18;

%% Compute spectral measure
M = size(PX,1);
TH2 = -pi:0.005:pi;                     % angles for spectral measure
epsilon = 0.1;                         % smoothing parameter
order = 2;                              % order of kernel

sig = PX(:,n_pendula+1)';
sig = sig - mean(sig);
g_coeffs = PX\transpose(sig); % cofficients of g in dictionary expansion

[~,xi] = riggedDMD(PX,PY,1/M,epsilon,[],[],'order',order,'g_coeffs',g_coeffs,'TH2',TH2);

figure % spectral measure plot
plot(TH2,xi,'linewidth',1)
xlim([-pi,pi])
ax=gca; ax.FontSize=18;
xlabel('$\theta$','interpreter','latex','fontsize',18)
% 
% 
% %% Correlation based method
% 
% MU = ErgodicMoments(sig,round(10/epsilon)); % computes correlations using ergodic thoerem - need to use Monte Carlo instead if using Monte Carlo sampling of trajectories
% mu1 = MomentMeas(MU,'filt','vand'); % given moments, computes measure with filter "filt" (see description of code)
% 
% figure
% plot(mu1,'linewidth',1)
% xlim([-pi,pi])
% xlabel('$\theta$','interpreter','latex','fontsize',18)
% ax = gca; ax.FontSize = 18;



