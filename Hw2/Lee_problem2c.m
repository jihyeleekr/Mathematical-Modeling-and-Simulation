function Lee_problem2c
% This function initializes arrays to store the positions of a random walker
% after 10, 40, and 160 steps. It then simulates random walks, where the
% walker takes steps based on specified probabilities. The walker's
% position is recorded after each specified number of steps, and histograms
% are plotted to visualize the distribution of positions.
%
% Parameters:
%   - samples: Number of random walk samples
%   - steps: Number of steps for histograms
% 
% Created
%   02/05/2024 by Jihye Lee


% Parameters
samples = 100000; % Number of random walk samples
steps = [10, 40, 160]; % Number of steps 

% Initialize 
positions_10 = zeros(1, samples); % Array for positions_10 steps
positions_40 = zeros(1, samples); % Array for positions_40 steps
positions_160 = zeros(1, samples); % Array for positions_160 steps

% Simulate random walks
for s = 1:samples %Time loop
    position = 0; % The starting position of an agent
    for step = 1:max(steps)
        % Determine step size based on probabilities
        if rand() <= 0.95 % with 95% of change the agent moves one step
            step_size = 1;
        else %1- 0.95 = 0.5
            step_size = 21; % with 5% of changes the agent moves 21 steps 
        end

        % Determine direction whether the agent moves up or down with
        % probability of 50% each
        direction = randi(2) * 2 - 3;

        % Update 
        position = position + direction * step_size; % Update the position of the agent based on the direction with step_size

        % After each 10,40, and 160 steps store the agents position into
        % corresponding arrays
        if step == 10
            positions_10(s) = position;
        end
        if step == 40
            positions_40(s) = position;
        end
        if step == 160
            positions_160(s) = position;
        end
    end
end

% Plot histograms
figure;

% Plot the histogram after 10 steps
subplot(3, 1, 1);
histogram(positions_10);
title('Histogram after 10 Steps');
xlabel('Position');
ylabel('Probability');

% Plot the histogram after 40 steps
subplot(3, 1, 2);
histogram(positions_40);
title('Histogram after 40 Steps');
xlabel('Position');
ylabel('Probability');

% Plot the histogram after 160 steps
subplot(3, 1, 3);
histogram(positions_160);
title('Histogram after 160 Steps');
xlabel('Position');
ylabel('Probability');

end
