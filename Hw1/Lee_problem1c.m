function problem1c

% Parameters
n = 40; % number of agents

% Initialization
x = double(1:n < 10); % The first five agents are standing up


% Computation
for j = 1:80 % time loop
    % Plotting
    fprintf('%d',x)
    fprintf('\n')

    %Update
    x = x([end 1:end-1]); % The wave moves to the right.
    
    pause % wait until press the key
end