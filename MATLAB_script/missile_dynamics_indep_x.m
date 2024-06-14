function State_dot = missile_dynamics_indep_x(x, State, u)
% This is missle dynamics with x as an independent variable.
% This doesn't solve mass (mass is already fixed as function of time)
% x: downrange (indep variable)
% State: state vectors, Here, XX = [t; h; v; gamma]
% u: control inputs. Here, u = alpha (AOA).

t = State(1); h = State(2); v = State(3); gamma = State(4);
alpha = u;

% units
% simulation time, t: sec
% horizontal position, x: m
% vertical positoin, h: m
% velocity, v: m/s
% flight path angle, gamma: rad
% mass, m: kg
% angle of attack, alpha: rad

% Universal constant
g = 9.81; % m/s2 (gravitational acc)

% Aerodynamics
[rho, Mach] = DensityMach(h, v);
CL_alpha = 8.5; % /rad, lift derivative
CD_0 = 0.30; % -, zero AOA drag coefficient
k = 0.02; % -, induced drag coefficient
S_ref = 0.026; % m2, reference area
CL = CL_alpha * alpha; % -, lift coefficient
CD = CD_0 + k * CL^2; % -, drag coefficient
q = 0.5 * rho * v^2; % N/m2, dynamic pressure
L = CL * q * S_ref; % N, lift
D = CD * q * S_ref; % N, drag

% Propulsion
tb = 7; % sec, burn time
mp = 60; % kg, propellant mass
m0 = 150; % kg, initial mass
if t <= tb
    T = 20000; % N, thrust before burnout
    mu = mp/tb; % kg/s, mass rate
    m = m0 - mu * t; % kg, mass at time t <= tb
else
    T = 0; % N, thrust after burnout
    mu = 0; % kg/s, mass rate
    m = m0 - mu * tb; % kg, mass at time t > tb
end

% EOM
dtdx = 1/(v*cos(gamma));
dhdx = tan(gamma);
dvdx = (T-D)/m/v/cos(gamma) - g/v * tan(gamma);
dgammadx = L/m/v^2/cos(gamma) - g/v^2;

State_dot = [dtdx; dhdx; dvdx; dgammadx];
end