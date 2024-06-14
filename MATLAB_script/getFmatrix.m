function F = getFmatrix(X_LG, U_LG, N, n, taus, func)
% calculate derivatives of the state
F = zeros(N, n);
for i = 1:1:N
    t_i = taus(i+1);
    tmp = func(t_i, X_LG(i, :)', U_LG(i, :));
    F(i, :) = tmp';
end

end