function Lee_problem4a
%TEMPLE_ABM_POPULATION_MIGRATE_MATE_AND_AGE
%   Agents (representing an animal species) move in a 2d
%   box. In each step an agent moves with a fixed speed
%   in its direction and slightly changes his direction.
%   Moreover, agents age and die when reaching their
%   maximum age. Whenever an integer cell is occupied by
%   exactly two agents, a new agent is born in the center
%   of that cell.
%
% 02/2016 by Benjamin Seibold
%            http://www.math.temple.edu/~seibold/
% 
% Modified by Jihye Lee
%               02/20/2024

% Parameters
N =200; % initial umber of agents
ns = 5000; % number of steps
ax = [0 10 0 10]; % problem domain
v = .02; % speed of agents (distance traveled per step) changed from .01 to .02
max_age = 100; % each agents lives 100 steps
p = 0.9;    % probability of producing new agent

% Initialization
cx = ax(1):ax(2); cy = ax(3):ax(4); % cell boundaries
X = [ax(1) + (ax(2) - ax(1)) / 2*rand(N,1),... % initial positions
    ax(3) + (ax(4) - ax(3))*rand(N,1)]; % of agents ((x, y) ∈ [0, 5] × [0, 10])
D = rand(N,1)*2*pi; % initial angles of direction of agent
A = randi(max_age,N,1); % age of agents
total_agents = zeros(ns,1); % total number of agetns


for j = 1:ns % loop over steps
    % Update positions and directions
    X(:,1) = X(:,1)+v*cos(D); X(:,2) = X(:,2)+v*sin(D); % move agents
    D = D+0.1*randn(size(D)); % change direction of motion
    A = A+1; % increase age of each agent by 1
    
    % Let agents bounce of walls
    ind = (X(:,1)<ax(1)&cos(D)<0)|... % who is hitting a wall
        (X(:,1)>ax(2)&cos(D)>0); % horizontally
    D(ind) = pi-D(ind); % reverse x-direction
    ind = (X(:,2)<ax(3)&sin(D)<0)|... % who is hitting a wall
        (X(:,2)>ax(4)&sin(D)>0); % vertically
    D(ind) = -D(ind); % reverse y-direction
    X(:,1) = min(max(X(:,1),ax(1)),ax(2)); % move agents outside of
    X(:,2) = min(max(X(:,2),ax(3)),ax(4)); % domain back onto boundary
    
    % Mating: if more than seven agent and less than ten agents in same cell, create new agent
    [occupied_cells,~,agents_in_cell] = unique(floor(X),'rows');
    agent_counts = accumarray(agents_in_cell,1);
    cells_with_multiple_agents = occupied_cells(...
        agent_counts == 2,:);

    offspring_reproduce = rand(size(cells_with_multiple_agents,1),1)<p; %probabilty of creating new offspring(1/4)
    N_new = sum(offspring_reproduce); % number of new agents
    X = [X;cells_with_multiple_agents(offspring_reproduce,:)+.5]; % new agents in cell centers
    D = [D;rand(N_new,1)*2*pi]; % random direction
    A = [A;zeros(N_new,1)]; % age of new agents is 0

    % Remove agents who are too old
    ind = A<= max_age; % agents who do not die
    X = X(ind,:); D = D(ind); A = A(ind); % keep those
    
    % Store the number of agents at j steps into an array
    total_agents(j) = length(A);
    
    % Plotting
    % clf
    % plot([1;1]*cx,ax([3,4])'*(cx*0+1),'k-',... % draw boundaries
    %     ax([1,2])'*(cy*0+1),[1;1]*cy,'k-') % of cells
    % hold on
    % ind = A<=max_age/3; % young agents
    % plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[0 .8 0])
    % ind = max_age/3<A&A<=max_age*2/3; % middle aged agents
    % plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[.8 .8 0])
    % ind = max_age*2/3<A; % old agents
    % plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[1 0 0])
    % hold off
    % axis equal xy, axis(ax)
    % xlabel('x'), ylabel('y') % axis labels
    % title(sprintf('Migrating, mating, and aging animals (%d agents)',...
    %     length(A)))
    % drawnow

    if ~numel(A), break, end % stop of no agents left
end
% To create a plot illustrating the number of agents against the step count as it increases.
% To create a plot illustrating the number of agents against the step count as it increases.
figure;
plot(total_agents)
ylim([0 max(total_agents)*1.1]) % set y-axis limit to start from 0
xlabel("Number of steps")
ylabel("Number of recoverd agents")
title('Number of Recovered Agents Over Time')
drawnow