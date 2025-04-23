function plot_drug_effects(stats, slicetab, perturbation, doses)
% This function plots the dose-dependent drug effects on measures of cell
% speed, proliferation and migration, derived from slice and single cell 
% statistics.
%
% @authors: Madeleine Skeppås
% @date: 140624

M = [];
pert_tab = stats.(perturbation);
for i=1:height(pert_tab)
    M(:,i) = table2array(pert_tab{i}(:,1));
end

figure;
% Normalize each row with the strongest per-metric effect
M=M./max(abs(M'))';
M = M([1 5 2 3 4],:);

h=heatmap(M);
h.YDisplayLabels = {'cell speed', 'adMAD', 'D', 'alpha', 'growth rate'};
h.XDisplayLabels=stats.Row;
colormap(redbluecmap)
clim([-1 1]*1.2)
title(['Estimates of dose-dependent ' char(perturbation) ' effect on metrics'])
fontsize('scale', 1.5)

tbl_ext = slicetab;
tbl_ext.speed = cell(height(tbl_ext),1);

for i=1:height(tbl_ext)
    % Add speed
    tbl_ext.speed{i} = calculate_cell_speed(tbl_ext.traxs{i}, tbl_ext.trays{i},tbl_ext.delta_t(i));
end

c = [0.1804, 0.3961, 0.5490;  % Deep Hokusai Blue  
     0.8784, 0.4745, 0.4039;  % Veronese Coral Red  
     0.3019, 0.6275, 0.8509;  % Soft Sky Blue  
     0.9255, 0.2275, 0.3196;  % Bold Scarlet Red  
     0.6078, 0.8039, 0.9333;  % Light Powder Blue  
     0.9529, 0.6235, 0.4812;  % Salmon Pink  
     0.5078, 0.0078, 0.2784;  % Dark Burgundy Red
     0,      0,      0;
];

cellines = unique(tbl_ext.HGCC);
tbl_ext.perturbation = categorical(tbl_ext.perturbation);
tbl_ext.perturbation = reordercats(tbl_ext.perturbation, string([{'control'} setdiff(tbl_ext.perturbation, {'control'})']));
tbl_ext.dose = categorical(tbl_ext.dose);

figure;
for i=1:length(cellines)
    subplot(2,round(length(cellines)/2),i)
    tbl = tbl_ext(logical(tbl_ext.HGCC == string(cellines{i})),:);
    tbl.dose = categorical(tbl.dose);
    tbl.dose = removecats(tbl.dose, string(setdiff(unique(tbl_ext.dose),unique(tbl.dose))));

    h = daboxplot(cell2mat(tbl.speed),'groups',tbl.dose,...
        'xtlabels', unique(tbl.dose),'colors',c,'whiskers',0,...
        'scatter',1,'scattersize',30,'scatteralpha',0.5,...
        'boxspacing',1, 'boxwidth', 2); 
    
    % ylim([0 15])
    ylabel('ROI-level avg cell speed (um/h)')
    xlabel("Dose (uM)")
    title(char(cellines{i}))
    ax=gca;
    ax.Box = 1;
    ax.LineWidth = 1;
end

sgtitle(['Dose-dependent effect of ' char(perturbation) ' treatment'])
fontsize('scale',1.7)

% Calculate significance
p_values=[];
for i=1:length(cellines)
    hgcc = cellines(i);

    tbl = tbl_ext(logical(tbl_ext.HGCC == string(hgcc)),:);

    [p,~,stats] = anova1(cell2mat(tbl.speed), tbl.dose, "off");
    % multcompare(stats)
    p_values = [p_values; p];
end

figure;
plot_admad(slicetab, perturbation, doses)

figure;
plot_msd(slicetab, perturbation, doses)

end
