% my_fmincon function test
clear all; close all; clc;

fun = @(x) 1+x(1)/(1+x(2)) - 3*x(1)*x(2) + x(2)*(1+x(1));
c = @(x) 0;
d = @(x) [x(1) - 1; x(2) - 2];
rho = 10000;
lam0 = 0;
mu0 = [0; 0];

x0 = [100; 100];
[x_opt, f_val] = my_fmincon(fun, x0, c, d, lam0, mu0, rho)

[x_opt, f_val] = fmincon(fun, x0, [], [], [], [], [], [], @mycon)

function [c, ceq] = mycon(x)
    c = [x(1) - 1; x(2) - 2];
    ceq =[];
end