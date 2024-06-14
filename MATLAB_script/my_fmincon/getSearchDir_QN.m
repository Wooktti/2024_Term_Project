function [x_k_1, p_k, H_k_1] = getSearchDir_QN(func, H_k, x_k)
    % H_k: estimate of inverse Hessian at x_k
    
    % check positive definiteness of H_k
    
    g_k = calcGradient(func, x_k); % calculate gradient at x_k

    % if isPosDef
    %     p_k = -H_k*g_k; % calculate search direction
    % else
    %     p_k = -g_k; % if H_k is not posdef
    % end
    p_k = -H_k*g_k;

    % update H_k => H_k_1

    % perform line search to determine alpha_k
    func_line = @(alpha) func(x_k + alpha*p_k); % function values along p_k
    alpha_MIN = 0;
    alpha_MAX = 10;
    eps = 1e-12;
    MAX_ITER = 10000;
    [alpha_k, flag] = goldenSectionSearch(func_line, alpha_MIN, alpha_MAX, eps, MAX_ITER);
    if flag == 0
        error("goldenSectionSearch: exceed MAX_ITER")
    end

    s_k = alpha_k * p_k; % x_k_1 - x_k
    x_k_1 = x_k + s_k; % next search point
    g_k_1 = calcGradient(func, x_k_1); % gradient at next search point
    y_k = g_k_1 - g_k; % difference of gradient
    
    % DFP algorithms:
    H_k_1 = H_k - (1 / (y_k' * H_k * y_k)) * H_k * (y_k * y_k') * H_k' + (1 / (s_k' * y_k)) * (s_k * s_k');

end