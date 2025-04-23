function plot_msd(tbl, perturbation, doses_in)

tbl.interpolated_admad = cell(height(tbl),1);
tbl.interpolated_sumgreen = cell(height(tbl),1);
tbl = tbl(ismember(tbl.perturbation, ["control", perturbation]),:);
cellines = unique(tbl.HGCC);

deltat_min = min(tbl.delta_t);
perts = unique(tbl.perturbation);

linestyles = {"-", "--", ":", "-.", "-", "-"};
markers = {"none", "none", "none", "none", "none", "none"};
leg1=[];

for j=1:length(perts)

    tab_pert = tbl(tbl.perturbation == string(perts{j}),:);

    if(nargin > 1 && ~(perts{j} == "control"))
        doses = doses_in;
    else
        doses = unique(tab_pert.dose);
    end

    if(perts{j} == "control")
        col = [0.7,0.7,0.7];
    else
        normVals = (doses - min(doses)) / (max(doses) - min(doses));
        col = [145/255, 40/255, 54/255];
        lightRed = 0.3 * col + 0.7 * [1,1,1];
        col = (1 - normVals.') .* lightRed + normVals.' * col; 
    end
 
    for dose=1:length(doses)

        dosetab = tab_pert(tab_pert.dose == doses(dose),:);
    
        % Loop through the cell lines
        for i=1:length(cellines)
            hgcc = cellines{i};
            
            % Create the subset of stacks that belong to the current cell line
            tab = dosetab(dosetab.HGCC == string(hgcc),:);
    
            if(~isempty(tab))
                max_frames = max(tab.frames);
                deltat_max = max(tab.delta_t);
            
                t_max = 0:deltat_min:(max_frames*deltat_max);
            
                setz = unique(tab.set);
            
                sd_celline=nan(length(t_max),sum(cellfun(@(x) size(x,2), tab.traxs)));
            
                for k=1:(length(setz))
                    zet = setz(k);
                    subtable = tab(tab.set == zet,:);
            
                    delta_t = subtable.delta_t(1);
            
                    max_interpolated_length = length(0:deltat_min:(delta_t*subtable.frames(1)));
                    
                    % Create empty arrays to store all tracks from all stacks, based on the
                    % maximum number of frames in a stack, and the total number of tracks
                    % from all stacks. Arrays follow the same structure as structures traxs
                    % and trays, but are meant to store all such structures together.
                    % Rows: frames
                    % Columns : tracks
                    all_traxs = nan(max(subtable.frames), sum(cellfun(@(x) size(x,2), subtable.traxs)));
                    all_trays = nan(max(subtable.frames), sum(cellfun(@(x) size(x,2), subtable.traxs)));
                    
                    % Handle the first stack separately
                    xs = subtable.traxs{1};
                    ys = subtable.trays{1};
                
                    % Append the first tracks
                    all_traxs(1:height(xs),1:width(xs)) = xs;
                    all_trays(1:height(xs),1:width(xs)) = ys;
                
                    % This is where the appended tracks end
                    tailend = width(xs);
                            
                    % Iterate over the remaining stacks, appending tracks in the same way
                    for kk=2:height(subtable)
                        xs = subtable.traxs{kk};
                        ys = subtable.trays{kk};
                
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
            
                    % Append to the cell line table collecting all standard deviation
                    % values
                    start_idx = min(find(isnan(nanmean(sd_celline))));
                    sd_celline(1:height(sd_interpolated_set), start_idx:start_idx + width(sd_interpolated_set)-1) = sd_interpolated_set;
                end
            
                % Filter away any tracks shorter than 4 frames
                idx = sum(~isnan(sd_celline)) > 3;
                sd_celline = sd_celline(:,idx);
        
                % Convert squared displacement to um
                sd_um = sd_celline .* 1.85;
            
                t=0:deltat_min:(deltat_min*(height(sd_um)-1));
        
                % Plot mean + standard deviation (bold line thus becomes MSD)
                subplot(2,round(length(cellines)/2),i)
                hold on
                p = stdshade(t,sd_um', 0.1, col(dose,:), [hgcc ' (n = ' num2str(length(unique(tab.exp))) ')' ' (ROIs = ' num2str(height(tab)) ')' ' Perturbation: ' perts{j} ' Dose: ' num2str(doses(dose))], linestyles{i}, markers{i},8);
                leg1 = [leg1 p];
                xlabel('t (h)')
                ylabel('MSD (µm)^2')
                title(hgcc)
                ax=gca;
                ax.XLim = [0 80];
                ax.Box = 1;
                ax.LineWidth = 1;
            end
        end
    end
end

% Add figure details and adjust axes
% legend(leg1, 'Location', 'southoutside');
sgtitle('Cell migration (MSD)')
fontsize('scale',1.5)

end
