% Simple simulation, under no control (alpha = 0)
clear all; clc; close all;

load("prev_sol_N65.mat") % load X_LG, U_LG, taus, tf

t0 = 0;
tf = 100; % sec
dt = 1; % sec
t = t0:dt:tf; % sec
N = length(t);

% initial values
x0 = 0; % m
h0 = 10000; % m
v0 = 300; % m/s
gamma0 = 0; % rad
m0 = 150; % kg

X = zeros(5, N);
X(:, 1) = [x0; h0; v0; gamma0; m0]; % apply initial condition
u = zeros(1, N);
for i = 1:1:N-1
    t_curr = t(i);
    X_curr = X(:, i);
    u_curr = interp1(taus(2:end), U_LG, t_curr, 'linear', 'extrap');
    u(i) = u_curr;
    % RK4
    k1 = missile_dynamics(t_curr, X_curr, u_curr);
    k2 = missile_dynamics(t_curr + dt/2, X_curr + k1 * dt/2, u_curr);
    k3 = missile_dynamics(t_curr + dt/2, X_curr + k2 * dt/2, u_curr);
    k4 = missile_dynamics(t_curr + dt, X_curr + k3 * dt, u_curr);
    X(:, i+1) = X_curr + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

figure
plot(X(1, :), X(2, :))
xlabel("horizontal")
ylabel("vertical")

figure
plot(t, X(3, :))
xlabel("time")
ylabel("velocity")

figure
plot(t, X(4, :) * 180/pi)
xlabel("time")
ylabel("gamma [deg]")


figure
plot(t, X(5, :))
xlabel("time")
ylabel("mass [kg]")

figure
plot(t, u * 180/pi)
xlabel("time [sec]")
ylabel("AoA [deg]")