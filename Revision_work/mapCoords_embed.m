function [mappedInfo_total, mappedInfo_mean] = mapCoords_embed(X, Y, coordsX, coordsY, info, startidx)
% MAPCOORDS maps coordinates in (X,Y) to row indices in coordsX/coordsY,
% and fetches feature vectors from info{t} matrices.
%
% Inputs:
%   X, Y       - [time x cells] matrices of coordinates
%   coordsX    - {time x 1} cell array, each is [N_t x 1] vector of x coords
%   coordsY    - {time x 1} cell array, each is [N_t x 1] vector of y coords
%   info       - {time x 1} cell array, each is [N_t x features] matrix
%
% Outputs:
%   mappedInfo_total - {time x cells} cell array; each cell contains a 1×features vector
%                      (e.g. the row from info{t} corresponding to X(t,c),Y(t,c))
%   mappedInfo_mean  - [time x cells] matrix of means across the feature vector

[nt, nc] = size(X);
rowIdx = nan(nt, nc);

mappedInfo_total = cell(nt, nc);
mappedInfo_mean  = nan(nt, nc);

for t = 1:nt
    t_delayed = t + startidx - 1;
    cx = coordsX{t_delayed};
    cy = coordsY{t_delayed};
    infoMat = info{t_delayed};  % [N_t x features] matrix

    for c = 1:nc
        xval = X(t,c); 
        yval = Y(t,c);

        % Find matching coordinate
        match = find(cx == xval & cy == yval, 1);

        if ~isempty(match)
            rowIdx(t,c) = match;
        elseif t > 1 && ~isnan(rowIdx(t-1,c)) && ~isnan(xval)
            % fallback to previous valid index
            rowIdx(t,c) = rowIdx(t-1,c);
        else
            rowIdx(t,c) = NaN;
        end

        % Map feature vector
        if ~isnan(rowIdx(t,c))
            try
                mappedInfo_total{t,c} = infoMat(rowIdx(t,c), :);
                mappedInfo_mean(t,c)  = mean(mappedInfo_total{t,c});
            catch 
                mappedInfo_total(t,c) = mappedInfo_total(t-1,c); 
                mappedInfo_mean(t,c) = mappedInfo_mean(t-1,c);
            end
        else
            mappedInfo_total{t,c} = nan(1, size(infoMat, 2));
            mappedInfo_mean(t,c)  = NaN;
        end
    end
end
end
