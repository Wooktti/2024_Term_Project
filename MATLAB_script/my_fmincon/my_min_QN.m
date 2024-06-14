function [x_opt, opt_val, x_k_history] = my_min_QN(func, x0)
% minimize func(x) (unconstrained) with Quasi-Newton method
% x0: initial guess of the solution


n = length(x0); % dimension of x
x_curr = x0; % initialize current x
x_k_history = [x_curr];

H_curr = eye(n); % initialize estimate of inverse of Hessian
while true
    [x_next, p_curr, H_next] = getSearchDir_QN(func, H_curr, x_curr);
    eps = norm(x_next - x_curr);
    if eps < 1e-6
        x_opt = x_curr;
        opt_val = func(x_opt);

        % check PD-ness Hessian at solution candidate
        G_sol = calcHessian(func, x_opt); % gradient at converged x
        % g_sol = calcGradient(func, x_opt) % gradient at converged x

        if isPosDef(G_sol)
            % Hessian at converged point is positive definite
            % fprintf("eps: ");
            % disp(eps)
            break % sufficiency checked, and local minima found
        else
            % may have better direction
            % need to find better direction
            % G_sol may have negative eigenvalue
            % choose eigenvector corresponding to neg eigenvalue.
            [V, D] = eigs(G_sol);
            d = eigs(G_sol);
            idx = find(d < 0);
            p_new = V(:, idx); % choose eigvec corresponding to neg eigval
            func_line = @(alpha) func(x_opt + alpha*p_new); % function values along p_k
            alpha_MIN = 0;
            alpha_MAX = 10;
            eps = 1e-12;
            MAX_ITER = 10000;
            [alpha_k, flag] = goldenSectionSearch(func_line, alpha_MIN, alpha_MAX, eps, MAX_ITER);
            if flag == 0
                error("goldenSectionSearch: exceed MAX_ITER")
            end

            x_next = x_opt + alpha_k * p_new; % next trial point
        end

    end
    % update x and H
    x_curr = x_next;
    H_curr = H_next;
    x_k_history = [x_k_history, x_next];
end

end