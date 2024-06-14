function res = L_aug(x, fun, c_eq, d_ineq, lam, mu, rho)
    c = c_eq(x);
    d = d_ineq(x);
    res = fun(x) + lam' * c + rho * c' * c + rho * max(d + mu/2/rho, 0)' * max(d + mu/2/rho, 0) - mu'*mu / 4 / rho;
end