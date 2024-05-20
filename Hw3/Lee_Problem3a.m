function Lee_Problem3a

% Lee_problem3a

% This MATLAB code simulates the hill-topping behavior of butterflies in a
% 2D landscape. The butterflies conduct random walks, choosing between moving
% to a neighboring cell with the highest elevation or a random neighboring
% cell. The landscape is continuously modified based on the butterfly's
% position, creating visible channels. The simulation visualizes the paths of
% multiple butterflies as they interact with and modify the landscape over time.
%
% Parameters such as the number of butterflies, random walk steps, and
% elevation functions are defined to control the simulation dynamics.
% 
%  Created 
%  02/019/2024 by Jihye Lee 
% 
%  Modification of 
%   temple_abm_bacteria_run_and_tumble_and_eat.m
%   temple_abm_butterfly_animation.m
%   temple_abm_butterfly_corridor_width.m
% 
%   02/2016 by Benjamin Seibold
%            http://www.math.temple.edu/~seibold/


% Parameters
N_butterflies = 250; % number of butterflies
ns = 150; % number of random walk steps
q = 0.4; % probability to move to neighbor cell with highest elevation
f = @(x,y) max(100-sqrt((x-30).^2+(y-30).^2),... % elevation function,
    50-sqrt((x-120).^2+(y-100).^2)); % consisting of two conical humps
ax = [0 150 0 150]; % domain for plotting
x0 = [85, 95]; % initial position of butterfly

px = ax(1):ax(2); % x-vector for concentration field
py = ax(3):ax(4); % y-vector for concentration field
[PX, PY] = meshgrid(px, py); % generate 2d position matrices

F = f(PX, PY); % elevation data
rx = [-1; -1; -1; 0; 0; 1; 1; 1]; % x-coordinate of neighbor cells of origin
ry = [-1; 0; 1; -1; 1; -1; 0; 1]; % y-coordinate of neighbor cells of origin

% Plot background image (elevation)
imagesc(px, py, F) % plot elevation as color map
axis equal xy, axis(ax)
hold on 


for i = 1:N_butterflies % loop over number of butterflies(the loop should end if the butterfly lands on the top)
    title(sprintf('Random walk over elevation map (butterfly %d/%d)', i, N_butterflies))
    X = x0; % set initial position
    for j = 1:ns
        if rand<q % with probability q
            nf = f(X(1), X(2)) - F(X(end,2)+1, X(end,1)+1); % elev. of neighboring cells
            [val,ind] = max(nf); % value and index with highest elevation

            if val < f(X(2)+1,X(1)+1) % if current cell is higher than all
                break % neighboring cell (hill top), stop random walk
            end
        
        else % otherwise, with probability 1-q
            ind = randi(8); % choose random index for neighbor cell
        end
        step = [rx(ind),ry(ind)]; % step to neighbor cell
        X(end+1,:) = X(end,:)+step; % take step and append matrix
        % Update concentration field for each butterfly separately
        L = 0.25 * exp(-0.05 * ((PX - X(end, 1)).^2 + ...
            (PY - X(end, 2)).^2)); % local influence of this butterfly
        F = max(F - L, 0); % reduce concentration by local influence

            
    end
    % Plot random walk path (just computed segment)
    plot(X(:,1),X(:,2),'k.-')
    drawnow
    pause(.05)

    clf
    imagesc(px, py, F) % plot elevation as color map
    hold on
    axis equal xy, axis(ax)
    xlabel('x'), ylabel('y') % axis labels
    plot(x0(1), x0(2), 'ro') % starting position
    drawnow
    
end

hold off
end
