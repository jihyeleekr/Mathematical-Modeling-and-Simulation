function Lee_problem2b
% Lee_problem2b
% Random walkers in a one-dimensional space follow a simple rule in each
% time step. Each agent can take either one step up or one step down,
% meaning its position is increased or decreased by 1. The top and bottom
% agents have an equal probability (0.5) of moving one step up or one step
% down. The middle agent's movement depends on its position relative to
% the top and bottom agents. If it is exactly in the middle, the
% probability of moving one step up or one step down is 0.5. Otherwise,
% the probability is 0.6.
%
% Importantly, any update step that would lead to agents colliding or
% crossing paths is rejected, resulting in dependent paths for the agents.
% A histogram is plotted to depict the position of the middle agent after
% 100 steps. The simulation is run for at least 10,000 random walk
% instances.

% Created 
%  02/05/2024 by Jihye Lee 

% Modification of 
%   02/05/2024 by Jihye Lee 
%              Lee_problem2a



% Parameters
n = 3; % number of agents
ns = 100; % number of steps
nsim = 10000; % number of random walk simulations

% Initialization
x = [-5,0,5]; % initial positions of agents
position = zeros(1, nsim);% middle agent's positions

% Computation
for k = 1:nsim % simulation loop
    x0 = x; % create copy of old state
    for j = 1:ns % time loop
        while 1 % start loop
            x_new = x; % create a new state to check for collision

            for i = 1:n % iterate over agents
                if i == 2 % middle agent
                    if abs(x(1) - x(3)) == 2 % check if the middle agent is exactly in between the other two agents
                        x_new(i) = x(i) + (randi(2, 1) * 2 - 3); % Then the probability of one step up or down is 0.5
                    else
                        nearest = sign(x(3) - x(1)); % determine the nearest agent, so that the middle agent move towrds the nearer agent
                        x_new(i) = x(i) + (rand() <= 0.6) * nearest; % move with a probability of 0.6 towards the nearer agent
                    end
                else % bottom and top agents
                    x_new(i) = x(i) + (randi(2, 1) * 2 - 3); % the probability of one step up or down is always 0.5
                end
            end

            if all(diff(x_new) > 0) %if no walker collision or crossing, 
                % diff([1,3,6]) = 2 , 3 2 is from 3-1 and 3 is from 6-3
                % all is the condition. So all x in diff(x) > 0 then berak
                x = x_new; % update the current state
                break %terminate loop
            end
        end
        x0 = x; % update current state
    end
    position(k) = x(2); % save the position of the middle agent
end

% Plot histogram of the position of the middle agent
clf
histogram(position);
xlabel('Position of Middle Agent after 100 steps');
ylabel('Probability');


