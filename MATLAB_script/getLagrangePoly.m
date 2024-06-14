function coeff = getLagrangePoly(taus, values)
% calculate lagrange polynomial that passes through (taus, values), and
% return coefficients of the polynomial.

% Input
% taus: interpolation points (x-value)
% values: function values at taus (y-value)

% Output
% coeff: coefficients of polynomial

N = length(taus); % number of interpolation points -> N-1th order polynomial
if length(values) ~= N
    error("dimension of interpolation points and values should match")
end

coeff = zeros(1, N);

for i = 1:1:N
    coeff_l_i = 1;

    for j = 1:1:N
        if j == i
            continue
        end
        
        coeff_l_i = conv(coeff_l_i, [1, -taus(j)]) / (taus(i) - taus(j));

    end

    coeff = coeff + values(i) * coeff_l_i;
    
end


end