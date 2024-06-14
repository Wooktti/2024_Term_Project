function res = calcLagrangePolyDiff(i, t, taus)
% This function calculate ldot_i(t) where interpolation points are given as
% taus.

% Note that i = 0 means the very first interpolation points, which is
% taus(1).

i = i + 1; % to fix the index shifting
N = length(taus); % number of interpolation points

res = 0;
for j=1:1:N
    if j == i
        continue
    end
    
    prod = 1;
    for k=1:1:N
        if (k == i) || (k == j)
            continue
        end
        prod = prod * (t-taus(k)) / (taus(i)-taus(k));
    end

    res = res + prod / (taus(i) - taus(j));
end

end