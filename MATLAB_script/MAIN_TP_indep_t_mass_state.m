% AE551 Term Project
% 20190357 Jaewook Shim
clear all;
clc;
close all;

%% Convert optimal control problem -> parameter optimization problem (LG pseudo-spectral method)

%%%%%%%%%% Problem Settings - You can change here %%%%%%%%%%
n = 4; % number of state variables
m = 1; % number of control inputs
N = 60; % number of collocation point.

% initial value of state (need to be 1 x n vector)
X0 = [0, 10000, 300, 0]; % x0 = 0 km, h0 = 10 km, v0 = 300 m/s, gamma0 = 0 deg

% load solution with smaller N for initial guess for x, u, and tf
load("prev_sol_N55_mass_known_indep_t.mat")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% guess for states using previous solution.
[col_pts, ] = lgwt(N, 0, tf);
col_pts = sort(col_pts);
x_guess = interp1(taus(2:end), X_LG(:, 1), col_pts, 'linear', 'extrap');
h_guess = interp1(taus(2:end), X_LG(:, 2), col_pts, 'linear', 'extrap');
v_guess = interp1(taus(2:end), X_LG(:, 3), col_pts, 'linear', 'extrap');
gamma_guess = interp1(taus(2:end), X_LG(:, 4), col_pts, 'linear', 'extrap');

mass = @(t) 150 - (60/7) * min(t, 7);
m_history = mass(col_pts);

% x_guess = linspace(10, 100000, N);
% h_guess = [linspace(10100, 50000, floor(N/2)), linspace(55000, 5010, N-floor(N/2))];
% v_guess = [linspace(500, 1500, floor(N/2)), linspace(1600, 600, N-floor(N/2))];
% gamma_guess = [linspace(5, 30, floor(N/2)), linspace(35, -40, N-floor(N/2))] * pi/180;
% m_guess = [linspace(150, 90, floor(N/5)), linspace(90, 90, N-floor(N/5))];
% X_LG_guess = [x_guess; h_guess; v_guess; gamma_guess; m_guess]'; % N x n, states

X_LG_guess = [x_guess, h_guess, v_guess, gamma_guess];

% guess for controls
% U_LG_guess = zeros(N, m); % N x m, controls at collocation locations
U_LG_guess = interp1(taus(2:end), U_LG, col_pts, 'linear', 'extrap');

% plot initial guess
plotInitialGuess(col_pts, x_guess, h_guess, v_guess, gamma_guess, m_history, U_LG_guess)

%% solve using fmincon
params0 = []; % parameter vector to be optimized; total Nn + Nm + 1 paramters.
for i = 1:1:N
    params0 = [params0; X_LG_guess(i, :)'; U_LG_guess(i, :)']; % reshape X_LG and U_LG to column vector
end
params0 = [params0; tf]; % reshape tf, params is (Nn + Nm + 1) size vector.

options = optimoptions(@fmincon,'MaxFunctionEvaluations',1000000, 'MaxIteration', 1000000);
[sol, fval] = fmincon(@(params) objFunc(params, X0, n, m, N), ...
                        params0, ...
                        [], ...
                        [], ...
                        [], ...
                        [], ...
                        [], ...
                        [], ...
                        @(params) pseudoSpectral_nonlincon(params, X0, n, m, N), ...
                        options)

% params should be reshaped to X = [X0; X_LG], U_LG, and tf
X_LG = zeros(N, n); % N x n, states
U_LG = zeros(N, m); % N x m, controls
for i = 1:1:N
    X_LG(i, :) = sol(1:n, 1)';
    U_LG(i, :) = sol(n+1:n+m, 1)';
    sol = sol(n+m+1:end);
end
tf = sol;
X = [X0; X_LG]; % (N+1) x n, states
[tau_collocation, weights] = lgwt(N, 0, tf); % weights: N weights are returned that correspond to collocation points.
taus = [0; sort(tau_collocation)]; % (N+1) points. tau(1:N) includes collocation points.

x_history = X(:, 1);
h_history = X(:, 2);
v_history = X(:, 3);
gamma_history = X(:, 4);
alpha_history = U_LG;

% calculate history from lagrange polynomial
% dt = 1;
% time_lagrange = 0:dt:tf;
% X_lagrange = getLagrangeApproxSol(X, taus, tf, n, N, dt);
% 
% x_history_lagrangePoly = X_lagrange(:, 1);
% h_history_lagrangePoly = X_lagrange(:, 2);
% v_history_lagrangePoly = X_lagrange(:, 3);
% gamma_history_lagrangePoly = X_lagrange(:, 4);
% m_history_lagrangePoly = X_lagrange(:, 5);

figure
plot(x_history/1000, h_history/1000);
hold on
% plot(x_history_lagrangePoly/1000, h_history_lagrangePoly/1000);
xlabel("x [km]")
ylabel("h [km]")

figure
plot(taus, x_history / 1000);
hold on
% plot(time_lagrange, x_history_lagrangePoly/1000)
xlabel("time [sec]")
ylabel("x [km]")

figure
plot(taus, h_history / 1000);
hold on
% plot(time_lagrange, h_history_lagrangePoly/1000)
xlabel("time [sec]")
ylabel("h [km]")

figure
plot(taus, v_history);
hold on
% plot(time_lagrange, v_history_lagrangePoly)
xlabel("time [sec]")
ylabel("v [m/s]")

figure
plot(taus, gamma_history * 180/pi);
hold on
% plot(time_lagrange, gamma_history_lagrangePoly * 180/pi)
xlabel("time [sec]")
ylabel("\gamma [deg]")

figure
plot(taus, [mass(0); m_history]);
hold on
% plot(time_lagrange, m_history_lagrangePoly)
xlabel("time [sec]")
ylabel("mass [kg]")

figure
plot(taus(2:end), alpha_history * 180/pi);
xlabel("time [sec]")
ylabel("\alpha [deg]")


%% Function Definitions
function [c, ceq] = pseudoSpectral_nonlincon(params, X0, n, m, N)
% params should be reshaped to X = [X0; X_LG], U_LG, and tf
X_LG = zeros(N, n); % N x n, states
U_LG = zeros(N, m); % N x m, controls
for i = 1:1:N
    X_LG(i, :) = params(1:n, 1)';
    U_LG(i, :) = params(n+1:n+m, 1)';
    params = params(n+m+1:end);
end
tf = params;
X = [X0; X_LG]; % (N+1) x n, states

% get collocation points, which is root of n-th order Legendre polynomial.
[tau_collocation, weights] = lgwt(N, 0, tf); % weights: N weights are returned that correspond to collocation points.
taus = [0; sort(tau_collocation)]; % (N+1) points. tau(1:N) includes collocation points.
D = getDmatrix(taus); % N x (N+1) matrix
F = getFmatrix(X_LG, U_LG, N, n, taus, @missile_dynamics_mass_known); % N x n matrix

X_final = X0 + weights' * F;
x_f = X_final(1);
h_f = X_final(2);

x_T = 100000; % m
h_T = 5000; % m

% constraint
ceq = reshape(F - D*X, 1, N*n);
ceq = [ceq, x_f - x_T, h_f - h_T];
c = [abs(U_LG) - 20 * pi/180];
end


function fval = objFunc(params, X0, n, m, N)
% params should be reshaped to X = [X0; X_LG], U_LG, and tf
X_LG = zeros(N, n); % N x n, states
U_LG = zeros(N, m); % N x m, controls
for i = 1:1:N
    X_LG(i, :) = params(1:n, 1)';
    U_LG(i, :) = params(n+1:n+m, 1)';
    params = params(n+m+1:end);
end
tf = params;
X = [X0; X_LG]; % (N+1) x n, states

% get collocation points, which is root of n-th order Legendre polynomial.
[tau_collocation, weights] = lgwt(N, 0, tf); % weights: N weights are returned that correspond to collocation points.
taus = [0; sort(tau_collocation)]; % (N+1) points. tau(1:N) includes collocation points.
F = getFmatrix(X_LG, U_LG, N, n, taus, @missile_dynamics_mass_known); % N x n matrix

X_final = X0 + weights' * F;
v_f = X_final(3);

fval = -v_f;
end
    

