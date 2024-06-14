function [x_opt, J_opt] = my_fmincon(objFun, x0, c_eq, d_ineq, lam0, mu0, rho)
% Created at 2024.06.13 by Jaewook Shim

% Description
% parameter optimization solver under nonlinear equality and inequality
% constraints. This solver uses Augmented Lagrangian Method (ALM).

% Input
% objFun: objective function to be minimized. objFun should return column
% vector.
% x0: initial guess for solution. This is column vector.
% c_eq: function describing equality constraint, i.e., c_eq(x) = 0. c_eq
% should return column vector.
% d_ineq: function describing inequality constraint, i.e., d_ineq(x) <= 0.
% d_ineq shoud return column vector.
% lam0: initial guess for lambda (EQ)
% mu0: initial guess for mu (INEQ)
% rho: coefficient of penalty term

% Output
% x_opt: optimal solution
% J_opt: optimal function value
MAX_ITER = 10000;

k = 0;
x_curr = x0;
lam_curr = lam0;
mu_curr = mu0;

while k < MAX_ITER

    % define augmented lagrangian
    % L_aug = @(x) objFun(x) + lam_curr' * c_eq(x) + rho * c_eq(x)' * c_eq(x) + rho * max(d_ineq(x) + mu_curr/2/rho, 0)' * max(d_ineq(x) + mu_curr/2/rho, 0) - mu_curr'*mu_curr / 4 / rho;
    % disp(L_aug(x_curr))
    % minimize L_aug
    x_next = fminsearch(@(x) L_aug(x, objFun, c_eq, d_ineq, lam_curr, mu_curr, rho), x_curr);

    % update multiplier
    lam_next = lam_curr + 2 * rho * c_eq(x_next)
    mu_next = max(mu_curr + 2 * rho * d_ineq(x_next), 0)

    % check convergence
    if norm(x_next - x_curr) < 1e-6
        x_opt = x_next;
        J_opt = objFun(x_opt);
        break
    end

    x_curr = x_next;
    lam_curr = lam_next;
    mu_curr = mu_next;
    k = k + 1;
end


end