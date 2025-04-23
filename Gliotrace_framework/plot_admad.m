function plot_admad(tbl, perturbation, doses_in)

tbl.interpolated_admad = cell(height(tbl),1);
tbl.interpolated_sumgreen = cell(height(tbl),1);
tbl = tbl(ismember(tbl.perturbation, ["control", perturbation]),:);
cellines = unique(tbl.HGCC);

deltat_min = min(tbl.delta_t);

perts = unique(tbl.perturbation);

leg1 = [];
linestyles = {"-", "--", ":", "-.", "-", "-"};
markers = {"none", "none", "none", "none", "none", "none"};

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
        col = [13/255, 59/255, 102/255];
        lightBlue = 0.3 * col + 0.7 * [1,1,1];
        col = (1 - normVals.') .* lightBlue + normVals.' * col; 
    end

    for dose=1:length(doses)

        dosetab = tab_pert(tab_pert.dose == doses(dose),:);

        for i=1:length(cellines)
            hold on
            hgcc = cellines{i};
    
            % Create the subset of stacks that belong to the current cell line
            tab = dosetab(dosetab.HGCC == string(hgcc),:);
            
            % Loop through stacks to interpolate values according to the smallest
            % deltat
            for k=1:height(tab)
                t2 = linspace(0,tab.delta_t(k)*(length(tab.adMAD{k})-1), length(tab.adMAD{k}));
                t1 = 0:deltat_min:t2(end);
                tab.interpolated_admad(k) = {interp1(t2,tab.adMAD{k},t1)};
        
            end
        
            % Create empty variable for storing adMAD values
            values = nan(height(tab),max(cellfun(@length, tab.interpolated_admad)));
        
            % Loop through the stacks of the current cell line, extract 
            % the values and remove outliers
            for m=1:height(tab)
                vals = tab.interpolated_admad{m}';
                [~, tfrm] = rmoutliers(vals, 'quartiles'); % Identify outliers
                vals(tfrm) = mean(vals); % Replace outliers with mean
                values(m,1:length(tab.interpolated_admad{m})) = vals*100;
            end
            
             % Cut if the time-resolved measures are longer than 70 frames
            % if(size(values,2)>70)
            %     values = values(:,1:70);
            % end
        
            t=0:deltat_min:(deltat_min*(width(values)-1));

            if(~isempty(tab))
                % Plot curves as mean + standard deviation
                subplot(2,round(length(cellines)/2),i)
                hold on
                p = stdshade(t,values, 0.1, col(dose,:), [hgcc ' (n = ' num2str(length(unique(tab.exp))) ')' ' (ROIs = ' num2str(height(tab)) ')' ' Perturbation: ' perts{j} ' Dose: ' num2str(doses(dose))], linestyles{i}, markers{i},8);
                leg1 = [leg1 p];
                xlabel('Δt (h)')
                ylabel('% pixels moving')
                title(hgcc)
                ax=gca;
                ax.XLim = [0 80];
                ax.Box = 1;
                ax.LineWidth = 1;
            end
        end
    end
end
sgtitle('Frame-wise pixel re-arrangement (adMAD)')
fontsize('scale',1.5)
% legend(leg1, Location ="bestoutside")

end
