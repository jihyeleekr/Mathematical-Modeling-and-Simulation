function Lee_Final_Project

% Agent Parameters
N = 500 ; % initial number of agents (1.576 million people lives in Philadelphia in 2021)
i = 0.019; % initial positivity rate
new = 0.01; % probability to add new infected agent every hour
v = 0.013; % speed of agents (distance traveled per step)

% Time Parameters
hours_in_day = 24; % consider each step an hour
ns = 115 * hours_in_day; % number of steps (hours)
max_infection = 6 * hours_in_day; % duration of infection before quarantine (removal)
quarantine_over = 14 * hours_in_day; % duration after which agent is no longer infectious, and is immune

% Domain Parameters
ax = [0 10 0 10]; % problem domain
cx = ax(1):ax(2); cy = ax(3):ax(4); % cell boundaries

% Initialization
X = [ax(1)+(ax(2)-ax(1))*rand(N,1), ax(3)+(ax(4)-ax(3))*rand(N,1)]; % initial positions of agents (spread throuhgout entire domain)
I = double(rand(N,1)<i); % infection status for each agent (0=negative, >0 is the amount of hours that they have been positive)
I = I .* floor(rand(N,1) * (max_infection / hours_in_day)); % give each agent that is initally positive a random duration of infection (under the amount before they quarantine)
D = rand(N,1)*2*pi; % initial angles of direction of agent
H = zeros(0,2); % array to hold the history of size of population that is infected, and the maximum size of the population
M = ones(N,1); % array holding mask type, 1,2,3,4,(5) for different mask types (initially holds 1's)

% Math to find base chance for agent in square shared by infected agent to become infected (will be modified by q_modifier)
r0 = 5.7; % amount of agents a single infected agent will infect
avg_pop_density = N / ((ax(2)-ax(1)) * (ax(4)-ax(3))); % average amount of people in cell to infect
avg_cells_visited = (max_infection * v) / 4; % average amount of cells to visit, divided by 4
p_base = r0 / (avg_pop_density * avg_cells_visited * max_infection); % max_infection is amount of opportunities to infect
p_modifier = [1, 0]; % no mask, cloth mask, surgical mask, respirator, immune-regardless of mask

% Initialize mask state variables for each agent
tempM = rand(N,1); % array of size N holding randomly generated numbers between 0 and 1
for j = 1:N % for all agents
    for k = 1:4 % for the 4 different mask types
        if tempM(j) <= p(k) % determine which agents have which masks based on the wearing percentages
            M(j) = k; 
            break
        end
    end
end


for j = 1:ns % loop over steps
    X(:,1) = X(:,1)+v*cos(D); X(:,2) = X(:,2)+v*sin(D); % move agents

    % Update positions and directions
    for k = 1:size(D)
        if I(k) < max_infection % if not being quarantined,
            D(k) = D(k) + 0.1*randn; % change direction of motion
        end
    end
    
    % Contain agents in domain: let agents "bounce" off walls
    ind = (X(:,1)<ax(1)&cos(D)<0)|... % who is hitting a boundary
        (X(:,1)>ax(2)&cos(D)>0); % horizontally
    D(ind) = pi-D(ind); % reverse x-direction
    ind = (X(:,2)<ax(3)&sin(D)<0)|... % who is hitting a boundary
        (X(:,2)>ax(4)&sin(D)>0); % vertically
    D(ind) = -D(ind); % reverse y-direction
    X(:,1) = min(max(X(:,1),ax(1)),ax(2)); % move agents outside of
    X(:,2) = min(max(X(:,2),ax(3)),ax(4)); % domain back onto boundary

    [occupied_cells,~,agents_in_cell] = unique(floor(X),'rows'); 
    cells_with_multiple_agents = occupied_cells(accumarray(agents_in_cell,1)>=2,:); % define cells with multiple agents

    for k = 1:size(cells_with_multiple_agents,1) % for every cell with adequate number of agents
        agents_in_this_cell = zeros(0,0);
        covid_present = 0;

        for l = 1:size(X,1) % check every agent to see if they are in the cell, and if they are, add them to a list and check to see if they have COVID
            if (X(l,1) > cells_with_multiple_agents(k,1) && X(l,1) - cells_with_multiple_agents(k,1) <= 1 && X(l,2) > cells_with_multiple_agents(k,2) && X(l,2) - cells_with_multiple_agents(k,2) <= 1) % to see if they are in a cell together
                agents_in_this_cell = [agents_in_this_cell;l];

                if I(l) >= 1 && I(l) <= max_infection 
                    covid_present = 1; % take note of if any agents have COVID (and have yet to quarantine)
                end
            end
        end

        if covid_present % if COVID is present in this cell
            for l = 1:size(agents_in_this_cell) % find the agents without COVID, and with a probability of the base COVID probability times the adjusted odds ratio 
                if ~I(agents_in_this_cell(l))>=1 && rand<(p_base * p_modifier(M(agents_in_this_cell(l)))) % of the specific type of mask they were using
                    I(agents_in_this_cell(l)) = 1; % give them COVID (and start the countdown before their quarantine)
                end
            end
        end
    end

    for k = 1:size(I) % check every agent
        if I(k) > quarantine_over % if quarantine is over 
            M(k) = 5; % set mask status to immune
            I(k) = 0;
        elseif I(k)>=1 % to see if they are infected
            I(k) = I(k) + 1; % if so, add a day to their time before quarantine
        end
    end

    ind = I <= max_infection; % Remove agents who have been infected for 5.6 days (quarantine)
    X = X(ind,:);D = D(ind); I = I(ind);M = M(ind); % "On average, symptoms showed up in the newly infected person about 5.6 days"
    
    
    current_size = size(X,1) - nnz(I >= max_infection);
    
    if rand<new % With certain probability,
        X = [X;ax(1)+(ax(2)-ax(1))*rand, ax(3)+(ax(4)-ax(3))*rand];D = [D;rand*2*pi];I = [I;1];M = [M;1]; %  add an agent that has COVID
    end

    H = [H;[j, nnz(I >= 1 & I <= max_infection), current_size]]; %add point to the history
    
    % Plotting
    % Subplot 1:
    clf
    subplot(1,2,1);
    plot([1;1]*cx,ax([3,4])'*(cx*0+1),'k-',... % draw boundaries
        ax([1,2])'*(cy*0+1),[1;1]*cy,'k-') % of cells
    hold on
    if nnz(I > max_infection) % display quarantining agents in the legend if there are any
        leg = legend(sprintf('Quarantining (%d)', nnz(I > max_infection)));
    else
        leg = legend("");
    end
    leg.Location = 'northeastoutside'; % set location of legend 
    
    % Plot points (agents)
    ind = ~I>=1 & M==1; % uninfected agents with no mask
    plot(X(ind,1),X(ind,2), '.','markersize',12,'color',[0.6 0.2 0.8], 'DisplayName', sprintf('Not Infected (%d)',nnz(ind)));    
    ind = M==5; % infected agents that have finished quarantine
    plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[0.9 0.4 0], 'DisplayName', sprintf('Immune (%d)',nnz(ind)))
    % quarantining agents are not displayed
    ind = I>=1 & I <= max_infection; % infected agents that have not yet quarantined
    plot(X(ind,1),X(ind,2),'.','markersize',12,'color',[0.9 0 0.1], 'DisplayName', sprintf('Infected (%d)',nnz(ind)))

    hold off
    axis equal xy, axis(ax)
    title('Spread of COVID on phili')
    subtitle(sprintf('(%d agents, %d days and %d hours)', current_size, (floor(j/24)), mod(j,24)))

    % Subplot 2
    subplot(1,2,2);
    hold on
    plot(H(:,1), H(:,2), '-', 'color',[0.9 0 0.1]); % graph to show how many agents are positive for COVID
    text(1,current_size + 5, sprintf("   Current population size (%d)", current_size),'color', '#555555') 
    yline(current_size, ':','color', '#555555'); % line to show how many agents are not quarantining
    axis([1 max(j,10) 0 max(H(:,3))]);
    xlabel('Hours'), ylabel('Infected People')
    
    drawnow
end
