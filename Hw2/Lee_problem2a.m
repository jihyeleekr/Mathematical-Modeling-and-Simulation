function Lee_problem2a
%Lee_problem2a
%   Random walkers in one space dimension. In each time step,
%   each walker steps either one step up or one step down
%   (i.e., its position is increased or decreased by 1).
%   In addition, any update step that would cause walkers to
%   collide or cross paths is rejected. Hence the paths of the
%   walkers are not independent.
%   Plotted is the probability distribution of the position of the
%   bottom agent after 100 random walk steps, using 10,000 random walk runs.

% Created 
%  02/05/2024 by Jihye Lee

% Modification of 
%   01/2016 by Benjamin Seibold
%            http://www.math.temple.edu/~seibold/
% 
% Modification Date:
%           02/05/2024 by Jihye Lee

% Parameters
n = 3; % number of walkers
ns = 100; % number of steps
nsim = 10000; % number of random walk simulations

% Initialization
x = (0:n-1)*5; % initial positions of walkers
position = zeros(1, nsim); %bottom agent's positions

% Computation
for k = 1:nsim % simulation loop
    x0 = x; % create copy of old state
    for j = 1:ns % time loop
        while 1 % start loop
            x = x0 + (randi(2, 1, n) * 2 - 3); % add random step to all walkers
            %randi(2,1,n) generate 3 random numbers between 1 and 2. 
            % So randi(2,1,n)* either 2 or 4.
            % Then -3 make is 1 or -1 so that 1 move forward and -1 move backward.

            if all(diff(x) > 0) % if no walker collision or crossing,
                % diff([1,3,6]) = 2 , 3 2 is from 3-1 and 3 is from 6-3
                % all is the condition. So all x in diff(x) > 0 then berak
                break % terminate loop
            end
        end
        x0 = x; % update current state
    end
    position(k) = x(1); % save the position of the bottom agent
end

% Plot histogram of the position of the bottom agent
clf
histogram(position);
xlabel('Position of Bottom Agent after 100 steps');
ylabel('Probability');
