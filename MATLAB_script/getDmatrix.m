function D = getDmatrix(taus)
N1 = length(taus); % this is N+1
N = N1-1;
D = zeros(N, N1); % D: N x (N+1) size. D_ki = diff of l_i(tau_k)
for i = 0:1:N % i: index for all grid points, 0 ~ N (N+1 points)
    for k = 1:1:N % k:  index for collocation points, 1 ~ N (N points)
        tau_k = taus(k+1);
        D(k, i+1) = calcLagrangePolyDiff(i, tau_k, taus); % D_ki = ldot_i(tau_k)
    end
end

end