function grad_j = calcGradient_j(func, x, j)
% return jth component of gradient of func at x
    grad = calcGradient(func, x);
    grad_j = grad(j);
end