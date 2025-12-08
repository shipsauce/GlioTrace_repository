function mappedInfo_total = mapCoords_embed(X, Y, coordsX, coordsY, info, startidx)
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

[nt, nc] = size(X);
rowIdx = nan(nt, nc);

mappedInfo_total = cell(nt, nc);

for t = 1:nt
    t_delayed = t + startidx - 1;
    cx = coordsX{t_delayed};
    cy = coordsY{t_delayed};
    try
    infoMat = info{t_delayed};  % [N_t x features] matrix
    catch
        infoMat = [];
    end

    for c = 1:nc
        xval = X(t,c); 
        yval = Y(t,c);

        % Find matching coordinate
        match = find(cx == xval & cy == yval, 1);

        if ~isempty(match)
            rowIdx(t,c) = match;
        else
            rowIdx(t,c) = NaN;
        end

        % Map feature vector
        if ~isnan(rowIdx(t,c))
            mappedInfo_total{t,c} = infoMat(rowIdx(t,c), :);
        elseif ~isnan(xval)
            mappedInfo_total(t,c) = mappedInfo_total(t-1,c); 
        else
            mappedInfo_total{t,c} = nan(1, 256);
        end
    end
end
end
