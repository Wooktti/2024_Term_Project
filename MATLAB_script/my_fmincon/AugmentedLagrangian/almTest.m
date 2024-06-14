% Simple problems from wikipedia and elsewhere for testing the Augmented Lagrangian optimization.
%
% Implementation tips:
% ** Nonlinear programming problems tend to converge more reliably when the initial point, x0, is
%   feasible (feasibility of the initial point, however, is not required by the algorithm).
% ** Another strategy that may be used to improve convergence on nonlinear problems it to initialize
%   the Lagrange multipliers using the Moore-Penrose pseudoinverse:
%       lambda0 = pinv([dce(x0), dci(x0)])*df(x0)
%   which is the default initialization when lambda0 is not provided by the user (but with slack
%   variables included).
% ** Nonlinear problems may also require tuning of the increment parameter, gma, to improve
%   convergence.

p = 4;

gma = 10;

if p == 1
    % maximize f(x,y) = x+y subject to x^2+y^2=1
    % solution: global maximum @ (x,y) = (sqrt(2)/2,sqrt(2)/2).
    problem = struct( ...
        'f', @(x)(-sum(x)), ...
        'df', @(x)([-1; -1]), ...
        'ce', @(x)(sum(x.^2)-1), ...
        'dce', @(x)(2*x) ...
    );

    x0 = randn(2,1);
    lambda0 = 0;
elseif p == 2
    % maximize f(x,y) = x^2*y subject to x^2+y^2-3=0
    % solution: global maxima @ (x,y) = (sqrt(2),1) & (-sqrt(2),1), local maximum @ (x,y) = (0,-sqrt(3))
    problem = struct( ...
        'f', @(x)(-x(1)^2*x(2)), ...
        'df', @(x)([-2*prod(x); -x(1)^2]), ...
        'ce', @(x)(sum(x.^2)-3), ...
        'dce', @(x)(2*x) ...
    );

    x0 = rand(2,1);
    lambda0 = 0;
elseif p == 3
    % minimize f(x,y) = x^2+2*y^2+2*x+8*y subject to -x-2*y+10<=0, x>=0, y>=0
    % solution: global minimum @ (x,y) = (4,3)
    problem = struct( ...
        'f', @(x)(x(1)^2+2*x(2)^2+2*x(1)+8*x(2)), ...
        'df', @(x)([2*x(1)+2; 4*x(2)+8]), ...
        'ci', @(x)(x(1)+2*x(2)-10), ...
        'dci', @(x)([1; 2]), ...
        'lb', [0; 0] ...
    );

    x0 = [5; 5];
    lambda0 = 0;
elseif p == 4
    % find the maximum entropy distribution of a six-sided die under constraints {x}>=0 and sum({x})=1
    % solution: global maximum @ {x} = 1/6
    % This problem has difficulty converging unless lambda0 and gma are set to a reasonable value.
    x0 = rand(6,1);
    x0 = x0/sum(x0);
    gma = 15;

    problem = struct( ...
        'f', @(x)(sum(x.*log(x+eps))), ...
        'df', @(x)(log(x+eps) + 1), ...
        'ce', @(x)(sum(x)-1), ...
        'dce', @(x)(ones(size(x0))), ...
        'lb', zeros(size(x0)) ...
    );

    lambda0 = []; % use default
elseif p == 5
    % maximize f(x,y,z) = x*y*z subject to x+y+z=1, x>=0, y>=0, z>=0
    % solution: (x,y,z) = (1/3,1/3,1/3)
    problem = struct( ...
        'f', @(x)(-prod(x)), ...
        'df', @(x)([-x(2)*x(3); -x(1)*x(3); -x(1)*x(2)]), ...
        'ce', @(x)(sum(x)-1), ...
        'dce', @(x)(ones(3,1)), ...
        'lb', zeros(3,1) ...
    );

    x0 = zeros(3,1);
    lambda0 = 0;
elseif p == 6
    % minimize f(x,y,z) = 4*y-2*z subject to 2*x-y-z=2, x^2+y^2=1
    % solution: (x,y,z) = (2/sqrt(13),-3/sqrt(13),-2+7/sqrt(13))
    problem = struct( ...
        'f', @(x)(4*x(2)-2*x(3)), ...
        'df', @(x)([0; 4; -2]), ...
        'ce', @(x)([2*x(1)-x(2)-x(3)-2, x(1)^2+x(2)^2-1]), ...
        'dce', @(x)([2, 2*x(1); -1, 2*x(2); -1, 0]) ...
    );

    x0 = randn(3,1);
    lambda0 = zeros(2,1);
elseif p == 7
    % minimize f(x,y) = (x-2)^2+2*(y-1)^2 subject to x+4*y<=3, x>=y
    % solution: (x,y) = (5/3,1/3)
    problem = struct( ...
        'f', @(x)([1, 2]*(x-[2; 1]).^2), ...
        'df', @(x)([2*(x(1)-2); 4*(x(2)-1)]), ...
        'ci', @(x)([-x(1)-4*x(2)+3, x(1)-x(2)]), ...
        'dci', @(x)([-1, 1; -4, -1]) ...
    );

    x0 = [1; 0];
    lambda0 = zeros(2,1);
    gma = 50;
elseif p == 8
    % minimize f(x,y,z) = (x-1)^2+2*(y+2)^2+3*(z+3)^2 subject to z-y-x-1=0, z-x^2>=0
    % solution: (x,y,z) = (0.12288,-1.1078,0.015100)
    problem = struct( ...
        'f', @(x)((x(1)-1)^2 + 2*(x(2)+2)^2 + 3*(x(3)+3)^2), ...
        'df', @(x)([2*x(1)-2; 4*x(2)+8; 6*x(3)+18]), ...
        'ce', @(x)(x(3)-x(2)-x(1)-1), ...
        'dce', @(x)([-1; -1; 1]), ...
        'ci', @(x)(x(3) - x(1)^2), ...
        'dci', @(x)([-2*x(1); 0; 1]) ...
    );

    x0 = [0; 1; 0];
    lambda0 = zeros(2,1);
    gma = 110;
else
    return;
end

options = struct( ...
    'gma', gma ...
);

[x, fval, lambda, kkt] = almSolve(problem, x0, lambda0, options);
