function problem1d
    % Parameters
    n = 40; % number of agents
    

    % Initialization
    x = double(1:n < 6); % The first five agents are standing up

    % Radomly set initial positions to keep stand up and keep sit down
    keep_stand = randi([1, n], 1, 1);
    keep_sit = randi([1,n],1,1);
    
    % Computation
    for j = 1:80 % time loop
        % Plotting
        fprintf('%d', x);
        fprintf('\n');
        
        % Update the position
        x(keep_stand) = 1; % keep_stand position always 1
        x(keep_sit) = 0; % keep_sit position always 0
        
         % Check the left neighbor of keep_stand
        if x(keep_stand-1) == 1
            x(keep_stand+1) = 1;
        else
            x(keep_stand+1) = 0;
        end
        
         % Check the left neighbor of keep_sit
        if x(keep_sit-1) ==1
            x(keep_sit+1) = 1;
        end

        % Update the array (move the wave to the right)
        x = x([end, 1:end-1]);

        pause % wait until press the key
    end
end



