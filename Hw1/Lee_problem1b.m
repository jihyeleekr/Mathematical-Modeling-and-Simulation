function problem1b

% Parameters
n = 40;

% Initialize
wave = zeros(1,n);
right = round(n/3); % wave moves to the right
left = round(n/15); % wave moves to the left

wave(right) = 1; % 1 means the wave moves to the right
wave(left) = 2; % 2 means the wave moves to the left

for j = 1:80 % time loop
    fprintf([repmat(' %1.0f ',1,numel(wave)) '\n'],wave); % print wave array

    wave(right) = 0; % Make the current standing person sit
    wave(left) = 0; % Make the current standing person sit

    right = right + 1; % The wave moves to the right
    left  = left - 1; % The wave moves to the left

    if right <= n % if new right, which is right + 1, less than n
        wave(right) = 1; % Then make a person who's sitting on the new right stand

    else % if new right is greater than n,
        right = 1; %  then make right starts from position 1
        wave(right) = 1;
    end

    if left >= 1 % if new left, which is left -1, greater than 1
        wave(left) = 2; %Then make a person who's sitting on the new left stand

    elseif left == 0 % if new right is less than 1
        left = n; % then make right starts from position n
        wave(left) = 2;
    end

    pause
end

