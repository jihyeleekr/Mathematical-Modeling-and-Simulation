function Lee_problem5c
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
n = 22; % number of vehicles
L = 230; % length of road
w = 2; % width of road (for plotting only)
lv = 4.5; % length of each vehicle
tf = 60*5; % final time
dt = 0.1; % time step (for integration)
np = 1; % number of time steps between two plotting events
V = @(d) 10*(tanh(d/2-2)+tanh(2))/(1+tanh(2)); % optimal velocity function
f = @(x) [x(n+1:2*n);... % ODE right hand side
    20*(x([n+2:2*n,n+1])-x(n+1:2*n))./([x(2:n);x(1)+L]-x(1:n)-lv).^2+...
    .5*(V([x(2:n);x(1)+L]-x(1:n)-lv)-x(n+1:2*n))];

% Initialization
q = linspace(0,L,n+1)'; q = q(1:end-1); % initial positions of vehicles
v = (1:n)'; % initial velocities of vehicles
x = [q;v]; % initial state vector
nt = ceil(tf/dt/np)*np; % number of time steps (muliple of np)
dt = tf/nt; % actual time step
phi = linspace(0,2*pi,101); cx = cos(phi); cy = sin(phi); % circle
r = L/(2*pi); % radius of road

% Computation
t = 0; % initial time
for j = 0:nt/np % time loop
    % Computation
    for i = 1:(np*(j>0)) % Do np compute steps (first time: do nothing)
        k1 = f(x); % first slope
        k2 = f(x + dt/2*k1); % second slope
        k3 = f(x + dt/2*k2); % thrid slope
        k4 = f(x+dt*k3); % foruth slope
        x = x + dt/6*(k1+2*k2+2*k3+k4); % forward Runge-Kutta 4 step
        t = t+dt; % adance time
    end
    % Plotting
    q = x(1:n); % vehicle positions
    qv = [q-lv,q]'*2*pi/L; % vehicle front and back position angle
    clf
    subplot(1,2,1) % plot vehicles on circular road
    plot((r+w)*cx,(r+w)*cy,'k-',(r-w)*cx,(r-w)*cy,'k-') % draw road
    hold on
    plot(r*cos(qv),r*sin(qv),'b-','linewidth',3)
    hold off
    axis equal tight
    title(sprintf('Car-following model at t=%0.0fs: positions',t))
    subplot(1,2,2) % plot velocities over positions on line
    plot(mod(q,L),x(n+1:2*n),'b.')
    axis([0 L 0 20])
    title(sprintf('Car-following model at t=%0.0fs: velocities',t))
    xlabel('vehicle front position [m]')
    ylabel('vehicle velocity [m/s]')
    drawnow
end
