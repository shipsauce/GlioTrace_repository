function [traX, traY, phenotypes] = connect_tracklets(traX, traY, phenotypes)
% This function iterates over the frames in the stack and tries to connect 
% track starts and stops that are within a given matching distance that grows
% iteratively. When a match is made between a stop and a start, the newest
% track (the start) is appended to the older track.
% 
% Input parameters:
% traX - frame-wise x-coordinates of tracks (rows = frames, cols = tracks)
% traY - frame-wise y-coordinates of tracks (rows = frames, cols = tracks)
% phenotypes - frame-wise properties of cells (1-5: morphology class probability
%                          scores, 6: morphology class, 7: TME interaction class)
%
% Output parameters:
% Same as input parameters but adjusted for connected tracks.
%
% @authors: Madeleine Skeppås, Sven Nelander
% @date: 14082024

    % combine tracklets into tracks
    maxdist = 7;
    maxdist2 = 15;
    
    % iterate over a range of distances for matching tracklets
    for radius=maxdist:maxdist2
        % iterate over the frames
        for i=3:size(traX,1)
            % find stops by identifying NaNs in the i:th row
            % the i-1:th row must be notnan in order for the stop to be valid
            track_stops=find((isnan(traX(i,:))).*(~isnan(traX(i-1,:))));
    
            % find starts by identifying notnans in the i:th row
            % the i-1:th row must be NaN in order for the start to be valid
            track_starts=find((~isnan(traX(i,:))).*isnan(traX(i-1,:)));
        
            % p is the x,y-coordinates of newly ended tracks
            p=[traX(i-1,track_stops)' traY(i-1,track_stops)'];
            
            % q is the x,y-coordinates of the newly begun tracks
            q=[traX(i,track_starts)' traY(i,track_starts)'];
            
            % calculate pairwise distances and pair ends of tracklets
            cost=pdist2(q,p);
            pairs=matchpairs(cost,radius);
            
            % iterate over the pairs and connect the newly begun tracks to the
            % end of the ended tracks, and erase the newest track
            for j=1:size(pairs,1)
                traX(i:end,track_stops(pairs(:,2)))=traX(i:end,track_starts(pairs(:,1)));
                traY(i:end,track_stops(pairs(:,2)))=traY(i:end,track_starts(pairs(:,1)));
                traX(i:end,track_starts(pairs(:,1)))=nan;
                traY(i:end,track_starts(pairs(:,1)))=nan;
                
                % apply the same operation to every phenotype
                for k=1:length(phenotypes)
                    data = phenotypes{k};
                    data(i:end,track_stops(pairs(:,2))) = data(i:end,track_starts(pairs(:,1)));
                    data(i:end,track_starts(pairs(:,1)))=nan;
                    phenotypes{k} = data;
                end
            end
        end
    
        % select only the tracks that have not been erased
        f=find(sum(isnan(traX))<size(traX,1));
        traX=traX(:,f);
        traY=traY(:,f);

        for p=1:length(phenotypes)
            data = phenotypes{p};
            data = data(:,f);
            phenotypes{p} = data;
        end

  
    end

end