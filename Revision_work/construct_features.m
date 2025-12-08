function [localDensity, microgliaFraction, polarization, vesselCloseness] = construct_features(x_coords, y_coords, tme_status, vesselMask, delta_t, radius)

[numTime, numCells] = size(x_coords); % Get number of timepoints and cells

% Preallocate arrays
localDensity = nan(numTime, numCells);
microgliaFraction = nan(numTime, numCells);
polarization = nan(numTime, numCells);
vesselCloseness = nan(numTime, numCells);

% Compute velocities between timepoints for each cell, using time-dependent
% deltat
dx = nan(numTime, numCells);
dy = nan(numTime, numCells);
for c = 1:numCells
    for t = 2:numTime
        if ~isnan(x_coords(t,c)) && ~isnan(x_coords(t-1,c))
            dx(t,c) = (x_coords(t,c) - x_coords(t-1,c)) / delta_t(t-1,c);
            dy(t,c) = (y_coords(t,c) - y_coords(t-1,c)) / delta_t(t-1,c);
        end
    end
end

% Compute average velocity of present cells at each timepoint
roi_dx = nan(numTime,1);
roi_dy = nan(numTime,1);
for t = 2:numTime
    present = ~isnan(dx(t,:));
    if any(present)
        roi_dx(t) = mean(dx(t,present));
        roi_dy(t) = mean(dy(t,present));
    end
end
roiSpeed = sqrt(roi_dx.^2 + roi_dy.^2); % Compute ROI speed in each timepoint

% Loop through timepoints, loop through present cells and check for
% neighbors, vessel closeness and polarization
for t = 1:numTime
    presentCells = find(~isnan(x_coords(t,:)));
    
    % Local density and vessel closeness
    for c = presentCells
        % Local density
        neighbors = 0;
        for c2 = presentCells
            if c2 ~= c
                dist = sqrt((x_coords(t,c)-x_coords(t,c2))^2 + (y_coords(t,c)-y_coords(t,c2))^2);
                if dist <= radius
                    neighbors = neighbors + 1;
                end
            end
        end
        localDensity(t,c) = neighbors; % Local density of each cell in each timepoint
        
        % For current present cell, check smallest distance to vessels
        [vx, vy] = find(vesselMask);
        if ~isempty(vx)
            dists = sqrt((x_coords(t,c) - vx).^2 + (y_coords(t,c) - vy).^2);
            vesselCloseness(t,c) = min(dists);
        end
    end
    
    % In each timepoint, check how many of the cells present are associated
    % with microglia relative to all present cells
    microgliaFraction(t,presentCells) = sum(tme_status(t,presentCells)==1) / length(presentCells);
    
    % For current timepoint, loop through present cells and copmute cosine
    % similarity of cells velocity vector and ROI mean velocity vector
    % 1 --> aligned with average ROI behavior
    % 0 --> orthogonal movement
    % -1 --> opposite direction
    if t >= 2
        for c = presentCells
            if ~isnan(dx(t,c)) && roiSpeed(t) > 0
                polarization(t,c) = (dx(t,c)*roi_dx(t) + dy(t,c)*roi_dy(t)) / ...
                                     (sqrt(dx(t,c)^2 + dy(t,c)^2) * roiSpeed(t));
            end
        end
    end
end

end
