function grad = calcGradient(func, x)
% calculate gradient of function at given point x
% this function uses central difference method
% x: n x 1 column vector
h = 1e-4;
n = length(x); 

grad = zeros(n, 1);
for i = 1:n
    ei = zeros(n, 1);
    ei(i, 1) = 1;

    grad(i, 1) = (func(x + h*ei) - func(x - h*ei))/2/h;
end

end