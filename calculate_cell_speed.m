function avgSpeed = calculate_cell_speed(x, y, deltaT)
% This function calculates the average speed across a time sequence
% for a single cell or a matrix of cells based on x and y 
% coordinates and time delay.
%
% Input parameters:
%   x       - Vector of x coordinates (1 x time)
%   y       - Vector of y coordinates (1 x time)
%   deltaT  - Time delay between consecutive time points (scalar)
%
% Output parameters:
%   avgSpeed - Average speed of the cell (scalar)
%
% @authors: Madeleine Skeppås
% @date: 16012025

% Number of cells
numCells = size(x,2);

% If x and y are matrices in the shape time x cells
if(numCells > 1)
    % Initialize an array to store speeds
    speeds = [];
    
    % Loop through each cell and calculate speed
    for cell = 1:numCells
        idx = ~isnan(x(:,cell));
        xs = x(idx,cell);
        ys = y(idx,cell);

        dx = diff(xs);
        dy = diff(ys);
        distances = sqrt(dx.^2 + dy.^2) * 1.85;
    
        % Compute speed as distance / time
        speed = distances / deltaT;
        
        speeds = [speeds mean(speed(~isnan(speed)))];
    end

    % Calculate the average speed
    avgSpeed = nanmean(speeds(~isnan(speeds)));

else
    % Compute the distance between consecutive time points
    dx = diff(x);
    dy = diff(y);
    distances = sqrt(dx.^2 + dy.^2) * 1.85;
    
    % Compute speed as distance / time
    speeds = distances / deltaT;
    
    % Calculate the average speed
    avgSpeed = mean(speeds);
end

end
