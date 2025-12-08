function [embds, hard_labs, x_coords, y_coords, tme_info, deltat_mat] = handle_missing_emissions_v2(embds, hard_labs, tme_info, x_coords, y_coords, deltat, labs)

% Loop through each cell and adjust deltat accordingly. Add everything
% back into matrices of the format time x cell. 

deltat_mat = nan(size(labs));

for cell=1:width(labs)
    startidx = min(find(~isnan(hard_labs(:,cell))));
    prop_lab = labs(:,cell);
    
    keep = ~prop_lab;
    keep = keep(startidx:end);
    deltat_cell = repmat(deltat, [1 sum(~isnan(hard_labs(:,cell)))-1])';
    dt_new = [];
    
    for i = 1:length(deltat_cell)
        if keep(i)
            last_valid = i;
            dt_new(end+1) = deltat_cell(i);
        else
            dt_new(last_valid) = dt_new(last_valid) + deltat_cell(i);
            dt_new(i) = 0;
        end
        
    end
    dt_new(dt_new==0) = NaN;
    deltat_mat(startidx:startidx + length(dt_new)-1, cell) = dt_new(1:end);
end

embds(logical(labs)) = {[]};
hard_labs(logical(labs)) = NaN;
tme_info(logical(labs)) = NaN;
x_coords(logical(labs)) = NaN;
y_coords(logical(labs)) = NaN;

% Sort out cells shorter than 3 tracks
tracklength = sum(~isnan(hard_labs),1);

idx = tracklength >= 3;
embds = embds(:,idx);
hard_labs = hard_labs(:,idx);
tme_info = tme_info(:,idx);
x_coords = x_coords(:,idx);
y_coords = y_coords(:,idx);
deltat_mat = deltat_mat(:,idx);

end