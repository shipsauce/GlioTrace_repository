function [traX, traY, x_hat_history, phenotypes, phenonames, Kn_history, startidx] = track_tumor_cells2(cellsx, cellsy, properties)
% This function takes cell coordinates from macro_track2 and associated
% properties from classify_tumor_cells and applies a Kalman filter to the
% cell coordinates to create single-cell trajectories.
% 
% Input parameters:
% cellsx - cell array with frame-wise vectors of x-coordinates
% cellsy - cell array with frame-wise vectors of y-coordinates
% properties - cell array with frame-wise tables of cell properties where
%               rows=cells and cols=features
%
% Output parameters:
% traX - matrix with time-resolved cell x-coordinates (row=frame, col=cell)
% traY - matrix with time-resolved cell y-coordinates (row=frame, col=cell)
% x_hat_history - matrix of historic state estimates (row=frame, col=cell)
% phenotypes - cell array of matrices of the time-resolved 
%                   properties of cells (row=frame, col=cell)
% phenonames - names of the properties in the phenotypes columns
% Kn_history - historic Kalman Gain estimates
% startidx - starting index of stack; the first frame where cells appear
%
% @authors: Madeleine Skeppås, Sven Nelander
% @date: 15082024

%% Define Kalman filter parameters
dt = 1; % Time step
nstates = 4;

F = [1, dt, 0, 0;
0, 1, 0, 0;
0, 0, 1, dt;
0, 0, 0, 1]; % State transition matrix

H = [1, 0, 0, 0;
0, 0, 1, 0]; % Measurement matrix

P = eye(4); % Initial state covariance

Q = 1 * eye(4); % Process noise covariance

R = 1 * eye(2); % Measurement noise covariance

% Find the first frame in the stack where cells first appear
try
    n_phenotypes=size(properties{1},2); % Count the number of object properties
    phenonames=properties{1}.Properties.VariableNames; % Store names of phenotypes
    startidx=1;
catch
    startidx = find(~cellfun(@isempty, properties),1,'first');
    n_phenotypes=size(properties{startidx},2); % Count the number of object properties
    phenonames=properties{startidx}.Properties.VariableNames; % Store names of phenotypes
end

x_hat = [cellsx{startidx}, zeros(length(cellsx{startidx}),1), cellsy{startidx}, zeros(length(cellsx{startidx}),1)]'; % Initial guess

% Initialize arrays to store the filter's estimated state and covariance
x_hat_history = {};
x_hat_history{1} = cellsx{startidx}';
x_hat_history{2} = zeros(length(cellsx{startidx}),1)';
x_hat_history{3} = cellsy{startidx}';
x_hat_history{4} = zeros(length(cellsx{startidx}),1)';
Kn_history = zeros(2, size(cellsx,2));
phenotypes={};
track_counter = dictionary(0,0);

%% Handle the first frame separately
% Extrapolate the first guess
x_hat_minus = F * x_hat;

P_minus = F * P * F' + Q;

traX=cellsx{startidx}';
traY=cellsy{startidx}';
current_pos = [cellsx{startidx} cellsy{startidx} table2array(properties{startidx})];

% Store phenotypes from the 1st frame in phenotypes object
for i=1:n_phenotypes
    phenotypes{i}=properties{startidx}.(i)';
end

%% Kalman filter loop

for j=startidx+1:length(cellsx)

    % Start with identifying cell pairs

    % Get cell positions and properties from the latest frame
    latest_props=zeros(size(traX,2),n_phenotypes);
    for ij=1:n_phenotypes
        latest_props(:,ij)=phenotypes{ij}(end,:)';
    end

    latest_pos=[traX(end,:)' traY(end,:)' latest_props];
    
    try
        current_pos=[cellsx{j} cellsy{j} table2array(properties{j})];
    catch
        current_pos = current_pos(:,1:2+n_phenotypes);
    end
    
    % Get the tracks of the latest position which are still monitored
    active_tracks=find(~isnan(latest_pos(:,1)));
    
    % Calculate the pairwise distances between predicted and measured positions
    cost=pdist2(single(x_hat_minus([1,3],:)'), single(current_pos(:,[1,2])));

    try
        % Elements in the cost matrix that are nan are set to 1000
        cost(find(isnan(cost)))=1000;
        % Pair the cells in the previous and current frame
        pairs=matchpairs(cost',18);
    catch
        1; 
    end

    % Extract indicies for paired cells in previous and current frame
    cells_paired_prev=pairs(:,2);
    cells_paired_curr=pairs(:,1);

    % Store the row numbers of current cells that weren't matched
    newtracks=setdiff(1:size(current_pos,1),cells_paired_curr);
    
    % Take care of old tracks that didn't get a match
    oldtracks = setdiff(active_tracks, cells_paired_prev);
    oldtracks_keep = [];

    for k=1:length(oldtracks)
        track = oldtracks(k);
        
        if(isKey(track_counter,track))
            track_counter(track) = track_counter(track) + 1;
        else
            track_counter(track) = 1;
        end
        
        if (track_counter(track) < 3)
            oldtracks_keep = [oldtracks_keep track];
        end
    end

    % Reset the count of tracks that were once again matched
    for mn=1:length(cells_paired_prev)
        track = cells_paired_prev(mn);
        
        if(isKey(track_counter,track))
            track_counter(track) = 0;
        end
    end

    % Now update the Kalman filter

    % Calculate the Kalman Gain Kn
    Kn = P_minus * H' / (H * P_minus * H' + R);
 
    % State update
    x_hat = nan(nstates, size(latest_pos,1));
    x_hat(:,cells_paired_prev) = x_hat_minus(:,cells_paired_prev) + Kn * (current_pos(cells_paired_curr,1:2)' - H * x_hat_minus(:,cells_paired_prev));
    x_hat(:,oldtracks_keep) = x_hat_minus(:,oldtracks_keep);

    % Covariance update (using the simplified and numerically unstable eq)
    P = (eye(4) - Kn * H) * P_minus;

    % Extrapolate the state
    x_hat_minus = F * x_hat;
    
    % Extrapolate the covariance
    P_minus = F * P * F' + Q;
    
    % Add state estimations to the information about the current frame, in
    % order to later store it properly.
    current_pos = [current_pos nan(size(current_pos,1), size(x_hat,1))];    
    
    for p=1:nstates
        current_pos(cells_paired_curr,n_phenotypes+2+p) = x_hat(p,cells_paired_prev)';            
    end

    % Extrapolate the positions of new tracks in order to intialize them
    newtrack_states = [current_pos(newtracks,1)'; zeros(1,length(newtracks)); current_pos(newtracks,2)'; zeros(1,length(newtracks))];
    x_hat_minus_newtracks = F * newtrack_states;
    x_hat_minus = [x_hat_minus x_hat_minus_newtracks];

    % Finally, we store the information about the matches

    % We select rows in the measurement table by choosing previous cells that
    % belonged to active tracks and then got paired
    % we fill these rows with data from current cells that got paired.
    measurement=nan(size(latest_pos,1),size(current_pos,2));
    % Measurement rows: cells, cols: properties
    measurement(cells_paired_prev,:) = current_pos(cells_paired_curr,:);
       
    % Measurement is vertically concatinated with the data from the current
    % unmatched cells, they are the new tracks.
    measurement=[measurement;current_pos(newtracks,:)];
    
    % Fill in the information about old tracks that will be kept.
    % Concatinate object information and predicted state from last frame.
    measurement(oldtracks_keep,:) = [latest_pos(oldtracks_keep,:) x_hat_minus(:,oldtracks_keep)'];

    % Make sure traX,traY have enough new columns to accomodate new tracks
    traX=[traX nan(size(traX,1),length(newtracks))];
    traY=[traY nan(size(traX,1),length(newtracks))];
    
    % Do the same thing for storage of historic estimates
    for kk=1:size(x_hat_history,2)
        x_hat_history{kk} = [x_hat_history{kk} nan(size(traX,1),length(newtracks))];
    end
    
    % Do the same thing for storage of phenotypes
    for k=1:n_phenotypes 
        phenotypes{k}=[phenotypes{k} nan(size(traX,1),length(newtracks))];
    end

    % Store the (x,y)-coordinates
    traX(end+1,:)=measurement(:,1)';
    traY(end+1,:)=measurement(:,2)';
    
    % Store the phenotypes
    for ii=1:n_phenotypes
        phenotypes{ii}(end+1,:)=measurement(:,2+ii)';
    end
    
    % Store the filters estimated state
    for m=1:nstates
        x_hat_history{m} = [x_hat_history{m}; measurement(:,n_phenotypes+2+m)'];
    end
    
    % Store the Kalman Gain for position and velocity
    Kn_history(1,j) = Kn(1);
    Kn_history(2,j) = Kn(2);

    %fprintf('Iteration: %d \n',j);
end

end