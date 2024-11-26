clc;
clear all;
close all;

addpath('scripts');

% Simulation parameters
n = 5; % Number of pendulums

L = 1.0;  % Total length of the chain of pendulums
m = 1.0;  % Total mass of the chain of pendulums

lengths = L/n * ones(n,1);  % Uniform lengths of pendulums
masses = m/n * ones(n,1);  % Uniform masses of pendulums
g = 9.81;  % Gravitational acceleration
t0 = 0; tf = 10; N_imp = 3000; N_ie = 10000; % Time range and steps

%Initial conditions
initial_theta = (7*pi/6) * ones(n, 1);
initial_p = zeros(n, 1);
z0 = [initial_theta ; initial_p];

XH = @(z) hamiltonianRHS(z, n, lengths, masses, g);
H = @(z) hamiltonianFunction(n, lengths, masses, g, z(1:n), z(n+1:end));

% Plot results

t_ie = linspace(t0,tf,N_ie+1);
t_imp = linspace(t0,tf,N_imp+1);

% Solve the system using the implicit midpoint method
[sol_imp, energy_imp] = implicit_midpoint_solver(XH, H, z0, t0, tf, N_imp);
[sol_ie, energy_ie] = implicit_euler_solver(XH, H, z0, t0, tf, N_ie);
[sol_eh, energy_eh] = euler_heun_solver(XH, H, z0, t0, tf, N_imp);
[sol_trpz, energy_trpz] = euler_heun_solver(XH, H, z0, t0, tf, N_imp);
[sol_ode45, t_ode45, energy_ode45] = ode45_solver(XH, H, z0, t0, tf, N_imp);

sols = {sol_imp,sol_ie,sol_eh,sol_trpz,sol_ode45};
energies = {energy_imp,energy_ie,energy_eh,energy_trpz,energy_ode45};
times = {t_imp,t_ie,t_imp,t_imp,t_ode45};
methods = {"Implicit Mid-Point","Implicit Euler","Euler Heun","Trapezoidal Method", "RK(5,4)"};
figure_title = sprintf('energy_behaviour_%d_pendula.pdf', n);
plotEnergies(energies, times, methods, figure_title)

%Create GIF for Implicit Mid-Point method
gif_filename = sprintf('chain_of_%d_pendula.gif', n);
% Generate the GIF
dt = t_imp(2)-t_imp(1);
generateGIF(sol_imp, lengths, n, dt, gif_filename, energy_imp);
