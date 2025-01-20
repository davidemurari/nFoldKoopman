addpath('scripts');

n_samples = 100;
n_pendula = 4;
n_steps = 11;
h = 0.05; %to generate a data point the method will do floor(h/0.01) substeps
plb = -2; %Lower bound for the interval where momenta are sampled
pub = 2; %Upper bound for the interval where momenta are sampled

[X,Y] = generateDataMontecarlo(n_samples, n_pendula, n_steps, h, plb, pub);
%X = sol_trpz(:,1:3000);
%Y = sol_trpz(:,2:3001);

%% MATRIX-BASED APPROACH

% %% Build matrix using kernel method
% [~,K,~,PX,PY] = kernel_ResDMD(X,Y,'type',"Gaussian",'N',500); % N = number of basis functions
% [~,LAM,W2] = eig(K,'vector');
% R = (sqrt(real(diag(W2'*L*W2)./diag(W2'*W2)-abs(LAM).^2))); % dual residual - warning, not L^2 residual

%% Build matrix using delay embedding
N = 10; % number of basis functions
M = 500; % number of data points
dim_sol = 2*n_pendula;
PX = zeros(M,N*dim_sol);
PY = zeros(M,N*dim_sol);
h = 1; % number of time steps for delay

for jj= 1:N
    PX(:,(jj-1)*dim_sol+1:jj*dim_sol) = transpose(X(:,(1:M)+(jj-1)*h));
    PY(:,(jj-1)*dim_sol+1:jj*dim_sol) = transpose(Y(:,(1:M)+(jj-1)*h));
end

K = PX\PY;
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

figure
loglog([0.001,1],[0.001,1],':k','linewidth',1)
hold on
loglog(sqrt(abs(abs(LAM).^2-1)),R,'k.','markersize',20)
xlabel('$\sqrt{1-|\lambda|^2}$','interpreter','latex','fontsize',18)
title('Dual Eigenpair Residuals','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%% Compute spectral measure
TH2 = -pi:0.005:pi;                     % angles for spectral measure
epsilon = 0.05;                         % smoothing parameter
order = 6;                              % order of kernel

sig = X(n_pendula+1,1:M);
sig = sig - mean(sig);
g_coeffs = PX\transpose(sig); % cofficients of g in dictionary expansion

[~,xi] = riggedDMD(PX,PY,1/(size(X,2)),epsilon,[],[],'order',order,'g_coeffs',g_coeffs,'TH2',TH2);

figure % spectral measure plot
plot(TH2,xi,'linewidth',1)
xlim([-pi,pi])
ax=gca; ax.FontSize=18;
xlabel('$\theta$','interpreter','latex','fontsize',18)


%% Correlation based method

MU = ErgodicMoments(sig,round(10/epsilon)); % computes correlations using ergodic thoerem - need to use Monte Carlo instead if using Monte Carlo sampling of trajectories
mu1 = MomentMeas(MU,'filt','vand'); % given moments, computes measure with filter "filt" (see description of code)

figure
plot(mu1,'linewidth',1)
xlim([-pi,pi])
xlabel('$\theta$','interpreter','latex','fontsize',18)
ax = gca; ax.FontSize = 18;



