% AE551 Term Project
% 20190357 Jaewook Shim
clear all;
clc;
close all;

%% Convert optimal control problem -> parameter optimization problem (LG pseudo-spectral method)

%%%%%%%%%% Problem Settings - You can change here %%%%%%%%%%
n = 4; % number of state variables
m = 1; % number of control inputs
N = 10; % number of collocation point.

% initial value of state (need to be 1 x n vector)
% indep variable: x from 0 to 100 km (fixed)
X0 = [0, 10000, 300, 0]; % t0 = 0 s, h0 = 10 km, v0 = 300 m/s, gamma0 = 0 deg

% load solution with smaller N for initial guess for x, u, and tf
load("prev_sol_N55_mass_known_indep_t.mat")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% guess for states using previous solution.
x_f = 100000; % m, final downrange
[col_pts, ] = lgwt(N, 0, x_f);
col_pts = sort(col_pts);
t_guess = interp1(X_LG(:, 1), taus(2:end) , col_pts, 'linear', 'extrap');
h_guess = interp1(X_LG(:, 1), X_LG(:, 2), col_pts, 'linear', 'extrap');
v_guess = interp1(X_LG(:, 1), X_LG(:, 3), col_pts, 'linear', 'extrap');
gamma_guess = interp1(X_LG(:, 1), X_LG(:, 4), col_pts, 'linear', 'extrap');

mass = @(t) 150 - (60/7) * min(t, 7);
m_history = mass(t_guess);

X_LG_guess = [t_guess, h_guess, v_guess, gamma_guess];

% guess for controls
% U_LG_guess = zeros(N, m); % N x m, controls at collocation locations
U_LG_guess = interp1(X_LG(:, 1), U_LG, col_pts, 'linear', 'extrap');

% plot initial guess
plotInitialGuess(t_guess, col_pts, h_guess, v_guess, gamma_guess, m_history, U_LG_guess)

%% solve using fmincon
params0 = []; % parameter vector to be optimized; total Nn + Nm paramters.
for i = 1:1:N
    params0 = [params0; X_LG_guess(i, :)'; U_LG_guess(i, :)']; % reshape X_LG and U_LG to column vector
end

lam0 = zeros(N*n + 1, 1);
mu0 = zeros(N, 1);
rho = 100;
[sol, fval] = my_fmincon(@(params) objFunc(params, X0, n, m, N), ...
                        params0, ...
                        @(params) c_eq(params, X0, n, m, N), ...
                        @(params) d_ineq(params, X0, n, m, N), ...
                        lam0, ...
                        mu0, ...
                        rho);

% params should be reshaped to X = [X0; X_LG], U_LG
X_LG = zeros(N, n); % N x n, states
U_LG = zeros(N, m); % N x m, controls
for i = 1:1:N
    X_LG(i, :) = sol(1:n, 1)';
    U_LG(i, :) = sol(n+1:n+m, 1)';
    sol = sol(n+m+1:end);
end
X = [X0; X_LG]; % (N+1) x n, states
[tau_collocation, weights] = lgwt(N, 0, x_f); % weights: N weights are returned that correspond to collocation points.
taus = [0; sort(tau_collocation)]; % (N+1) points. tau(1:N) includes collocation points.

x_history = taus;
t_history = X(:, 1);
h_history = X(:, 2);
v_history = X(:, 3);
gamma_history = X(:, 4);
alpha_history = U_LG;

figure
plot(x_history/1000, h_history/1000);
hold on
% plot(x_history_lagrangePoly/1000, h_history_lagrangePoly/1000);
xlabel("x [km]")
ylabel("h [km]")

figure
plot(t_history, x_history / 1000);
hold on
% plot(time_lagrange, x_history_lagrangePoly/1000)
xlabel("time [sec]")
ylabel("x [km]")

figure
plot(t_history, h_history / 1000);
hold on
% plot(time_lagrange, h_history_lagrangePoly/1000)
xlabel("time [sec]")
ylabel("h [km]")

figure
plot(t_history, v_history);
hold on
% plot(time_lagrange, v_history_lagrangePoly)
xlabel("time [sec]")
ylabel("v [m/s]")

figure
plot(t_history, gamma_history * 180/pi);
hold on
% plot(time_lagrange, gamma_history_lagrangePoly * 180/pi)
xlabel("time [sec]")
ylabel("\gamma [deg]")

figure
plot(t_history, [mass(0); m_history]);
hold on
% plot(time_lagrange, m_history_lagrangePoly)
xlabel("time [sec]")
ylabel("mass [kg]")

figure
plot(t_history(2:end), alpha_history * 180/pi);
xlabel("time [sec]")
ylabel("\alpha [deg]")


%% Function Definitions
function c = c_eq(params, X0, n, m, N)
% params should be reshaped to X = [X0; X_LG], U_LG, and tf
X_LG = zeros(N, n); % N x n, states
U_LG = zeros(N, m); % N x m, controls
for i = 1:1:N
    X_LG(i, :) = params(1:n, 1)';
    U_LG(i, :) = params(n+1:n+m, 1)';
    params = params(n+m+1:end);
end
X = [X0; X_LG]; % (N+1) x n, states

x_f = 100000; % m, final downrange (fixed)

% get collocation points, which is root of n-th order Legendre polynomial.
[tau_collocation, weights] = lgwt(N, 0, x_f); % weights: N weights are returned that correspond to collocation points.
taus = [0; sort(tau_collocation)]; % (N+1) points. tau(1:N) includes collocation points.
D = getDmatrix(taus); % N x (N+1) matrix
F = getFmatrix(X_LG, U_LG, N, n, taus, @missile_dynamics_indep_x); % N x n matrix

X_final = X0 + weights' * F;
h_f = X_final(2);

h_T = 5000; % m

% constraint
c = reshape(F - D*X, N*n, 1); % collocation constraint
c = [c; h_f - h_T]; % final alt: 5 km
end

function d = d_ineq(params, X0, n, m, N)
% params should be reshaped to X = [X0; X_LG], U_LG, and tf
X_LG = zeros(N, n); % N x n, states
U_LG = zeros(N, m); % N x m, controls
for i = 1:1:N
    X_LG(i, :) = params(1:n, 1)';
    U_LG(i, :) = params(n+1:n+m, 1)';
    params = params(n+m+1:end);
end
d = abs(U_LG) - 20 * pi/180;
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
X = [X0; X_LG]; % (N+1) x n, states

x_f = 100000; % m, final downrange

% get collocation points, which is root of n-th order Legendre polynomial.
[tau_collocation, weights] = lgwt(N, 0, x_f); % weights: N weights are returned that correspond to collocation points.
taus = [0; sort(tau_collocation)]; % (N+1) points. tau(1:N) includes collocation points.
F = getFmatrix(X_LG, U_LG, N, n, taus, @missile_dynamics_indep_x); % N x n matrix

X_final = X0 + weights' * F;
v_f = X_final(3);

fval = -v_f;
end
    

