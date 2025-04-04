function tbl_in = fit_msd(tbl_in)
cellines = unique(tbl_in.HGCC);

% Calculate MSD curve parameters for every combination of cell line, perturbation & dose, set

% Iterate over the cellines
for i=1:length(cellines)
    hgcc = cellines{i};
    subtable = tbl_in(tbl_in.HGCC == string(hgcc),:);

    perturbations = unique(subtable.perturbation);
    
    % Iterate over the different perturbations tested in celline
    for j=1:length(perturbations)
        pert = perturbations{j};
        
        subtable_1 = subtable(subtable.perturbation == string(pert),:);

        doses = unique(subtable_1.dose);
        
        % Iterate over the doses of current perturbation
        for jj=1:length(doses)
            
            dose = doses(jj);
            
            subtable_2 = subtable_1(subtable_1.dose == dose,:);
            
            idx = logical((tbl_in.HGCC == string(hgcc)) .* (tbl_in.perturbation == string(pert)) .* (tbl_in.dose == dose));
            
            max_frames = max(subtable_2.frames);
            deltat_max = max(subtable_2.delta_t);
            deltat_min = min(subtable_2.delta_t);
            
            t_max = 0:deltat_min:(max_frames*deltat_max);
            
            setz = unique(subtable_2.set);
            
            sd_subset=nan(length(t_max),sum(cellfun(@(x) size(x,2), subtable_2.traxs)));
            
            % Iterate over all the sets and interpolate the values from
            % each set separately, since each set can have distinct delta t
            
            for k=1:(length(setz))
                zet = setz(k);
                subtable_3 = subtable_2(subtable_2.set == zet,:);
                
                delta_t = subtable_3.delta_t(1);
                
                max_interpolated_length = length(0:deltat_min:(delta_t*subtable_3.frames(1)));
                
                % Create empty arrays to store all tracks from all stacks, based on the
                % maximum number of frames in a stack, and the total number of tracks
                % from all stacks. Arrays follow the same structure as structures traxs
                % and trays, but are meant to store all such structures together.
                % Rows: frames
                % Columns : tracks
                all_traxs = nan(max(subtable_3.frames), sum(cellfun(@(x) size(x,2), subtable_3.traxs)));
                all_trays = nan(max(subtable_3.frames), sum(cellfun(@(x) size(x,2), subtable_3.traxs)));
                
                % Handle the first stack separately
                xs = subtable_3.traxs{1};
                ys = subtable_3.trays{1};
                
                % Append the first tracks
                all_traxs(1:height(xs),1:width(xs)) = xs;
                all_trays(1:height(xs),1:width(xs)) = ys;
                
                % This is where the appended tracks end
                tailend = width(xs);
                
                % Iterate over the remaining stacks, appending tracks in the same way
                for kk=2:height(subtable_3)
                    xs = subtable_3.traxs{kk};
                    ys = subtable_3.trays{kk};
                    
                    all_traxs(1:height(xs),tailend+1:tailend+width(xs)) = xs;
                    all_trays(1:height(xs),tailend+1:tailend+width(xs)) = ys;
                    
                    tailend = tailend + width(xs);
                end
                
                % Create empty array for squared displacement of tracks
                sd = nan(size(all_traxs));
                
                % Iterate over all tracks (column by column)
                for m=1:width(all_traxs)
                    reference_pos_x=all_traxs(min(find(~isnan(all_traxs(:,m)))),m); % Find starting pos x of track
                    reference_pos_y=all_trays(min(find(~isnan(all_trays(:,m)))),m); % Find starting pos y of track
                    sd(1:length(find(~isnan(all_traxs(:,m)))),m) = abs(all_traxs(find(~isnan(all_traxs(:,m))),m) - reference_pos_x).^2 + ...
                    abs(all_trays(find(~isnan(all_trays(:,m))),m) - reference_pos_y).^2; % Calculate how far from the starting point each coordinate in track is
                end
                
                sd_interpolated_set = nan(max_interpolated_length, width(sd));
                
                % Loop through sd values to interpolate according to the smallest deltat
                for m=1:width(sd)
                    track = sd(~isnan(sd(:,m)),m);
                    tracklength = sum(~isnan(track));
                    t2 = linspace(0,delta_t*(tracklength-1), tracklength);
                    t1 = 0:deltat_min:t2(end);
                    try
                        interpolated_values = interp1(t2,track,t1);
                        sd_interpolated_set(1:length(interpolated_values),m) = interpolated_values;
                    catch
                        sd_interpolated_set(1:tracklength,m) = track;
                    end
                end
                
                % Append to the cell line table collecting all squared displacement
                % values
                start_idx = min(find(isnan(nanmean(sd_subset))));
                sd_subset(1:height(sd_interpolated_set), start_idx:start_idx + width(sd_interpolated_set)-1) = sd_interpolated_set;
            end
            
            % Convert squared displacement to um
            sd_um = sd_subset .* 1.85;

            t=0:deltat_min:(deltat_min*(height(sd_um)-1));

            % Calculate MSD
            msd = nanmean(sd_um,2);
         
            % Fit a straight line to get the diffusion 
            % parameters for the first time points
            log_y = log(msd(3:15));
            log_x = log(t([3:15]));
            p = polyfit(log_x, log_y, 1);
            log_y_fit = polyval(p, log_x);
            
            % Add information on D and alpha to the relevant ROIs for
            % current combination of celline, perturbation and dose
            tbl_in.alpha(idx) = {p(1)};
            tbl_in.D(idx) = {p(2)};
        
        end
        
    end

    i

end
end