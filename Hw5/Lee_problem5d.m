function Lee_problem5d
%TEMPLE_ABM_TRAFFIC_CAR_FOLLOWING
%   Model for vehicles driving on a single-lane road.
%   Each vehicle is following its immediate leader, and
%   accelerated and decelerates based on the headway, its
%   own velocity, and its leader's velocity. The vehicles
%   drive on a circular road. In this model, each vehicle
%   has two objectives: a) equilibrating its own velocity
%   to its leader's velocity; and b) having its own
%   velocity approach an optimal velocity, which is based
%   on its headway.
%
% 03/2016 by Benjamin Seibold
%            http://www.math.temple.edu/~seibold/
%Modified on
% 03/2024 by Jihye Lee

% Parameters
C = 0:0.05:2; % Array holding the coefficient in front of optimal velcoty term
n = 22; % number of vehicles
L = 230; % length of road
lv = 4.5; % length of each vehicle
tf = 60*5; % final time
dt = 1e-2; % time step (for integration)
np = 10; % number of time steps between two plotting events
V = @(d) 10*(tanh(d/2-2)+tanh(2))/(1+tanh(2)); % optimal velocity function
f = @(x,c) [x(n+1:2*n);... % ODE right hand side
    20*(x([n+2:2*n,n+1])-x(n+1:2*n))./([x(2:n);x(1)+L]-x(1:n)-lv).^2+...
    c*(V([x(2:n);x(1)+L]-x(1:n)-lv)-x(n+1:2*n))];

% Initialization
q = linspace(0,L,n+1)'; q = q(1:end-1); % initial positions of vehicles
v = (1:n)'; % initial velocities of vehicles
x = [q;v]; % initial state vector

nt = ceil(tf/dt/np)*np; % number of time steps (muliple of np)
dt = tf/nt; % actual time step

UMIN = zeros(n,1); % column vector to store minimum veloctiy at final time tf
UMAX = zeros(n,1); % column vector to store maximum veclotiy at final time tf


for s = 1:length(C)
    c = C(s); % Set the coefficient c to be desired value from array C

    % Computation
    t = 0; % initial time
    for j = 0:nt/np % time loop
        % Computation
        for i = 1:(np*(j>0)) % Do np compute steps (first time: do nothing)
            k1 = f(x,c); % first slope
            k2 = f(x + dt/2*k1,c); % second slope
            k3 = f(x + dt/2*k2,c); % thrid slope
            k4 = f(x+dt*k3,c); % foruth slope
            x = x + dt/6*(k1+2*k2+2*k3+k4); % forward Runge-Kutta 4 step
            t = t+dt; % adance time
        end
    end

    % Record umin and umax vlaues at final time into respective arrays
    UMIN(s) = min(x(n+1:2*n)); % min velocity in state vector
    UMAX(s) = max(x(n+1:2*n)); % max velocity in state vector
end

% Plotting a values vs UMIN and UMAX
clf
axis([0 L 0 20])
plot(C, UMIN, 'b-',C,UMAX,'g-');
legend('umin','umax')
title(sprintf('Investigating How Traffic Wave Characteristics Depend on the Strength of the Optimal Velocity Term'))
xlabel('Coefficient of velocity term [c]')
ylabel('Vechicle velocty [m/s]')

hold on
