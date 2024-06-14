function X_lagrange = getLagrangeApproxSol(X, taus, tf, n, N, dt)
% return smooth curves for solution, which is indeed lagrange polynomial
% X: [X0; X_LG], state values at interpolation points
% U_LG: control values at collocation points
% taus: interpolation points (tau0 + collocation points)
% tf: final time
% N: number of collocation points

time = 0:dt:tf; % sec
time = time';
N_time = length(time);

X_lagrange = zeros(N_time, n);
TimeMat = zeros(N_time, N+1);
for i=1:1:(N+1)
    TimeMat(:, i) = time.^(N-i+1);
end

% Calculate coefficient of lagrange basis polynomial
for j=1:1:n
    coeff = getLagrangePoly(taus, X(:, j))
    X_lagrange(:, j) = TimeMat * coeff';
end

end