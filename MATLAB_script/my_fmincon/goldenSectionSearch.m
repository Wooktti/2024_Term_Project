function [x_opt, flag] = goldenSectionSearch(func, a, b, eps, MAX_ITER)
% find x_opt = argmin func in [a, b]
% func is assumed to be unimodal
% this function uses 'golden section search'
% if flag is 1, the convergence occured within MAX_ITER.
% if flag is 0, the convergence didn't occur within MAX_ITER.

tau = (sqrt(5)-1)/2; % 0.6180

cnt = 0;
while ((b - a) > eps) && (cnt < MAX_ITER) % step 4. check uncertainty interval is reduced sufficiently small
    
    % step 1. evaluate f(x) at {0, 1-tau, tau, 1}
    f1 = func(a);
    f2 = func(a + (b-a)*(1-tau));
    f3 = func(a + (b-a)*tau);
    f4 = func(b);
    
    % step 2. determine which interval contains x_opt
    % step 3. select next interval and perform search.
    % A = [a, a + (b-a)tau]
    % B = [a + (b-a)(1-tau), b]
    if f2 <= f3
        % x_opt is in interval A
        b = a + (b-a)*tau; % update b to search in interval A
    else
        % x_opt is in interval B
        a = a + (b-a)*(1-tau); % update a to search in interval B
    end
    
    cnt = cnt + 1;
end

x_opt = (a + b)/2;

if cnt == MAX_ITER
    flag = 0;
else
    flag = 1;
end

end