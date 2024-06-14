function G = calcHessian(func, x)
% calculate hessian of scalar function 'func' at x
n = length(x);

G = zeros(n, n);
for j = 1:n
    gi_x = calcGradient(@(xx) calcGradient_j(func, xx, j), x);
    G(j, :) = gi_x';
end

end