function plot_drug_effects(statstab, perturbation)
% This function plots the dose-dependent drug effects on measures of cell
% speed and proliferation, derived from slice and single cell statistics.
%
% @authors: Madeleine Skeppås
% @date: 140624

figure;
stats = statstab.(char(perturbation));
M=table2array(stats(:,1));
M=M./max(abs(M'))';
M = M([1 5 2 3 4],:);

h=heatmap(M);
h.YDisplayLabels = {'cell speed', 'adMAD', 'D', 'alpha', 'growth rate'};
% h.XDisplayLabels=cellines(I);
colormap(redbluecmap)
clim([-1 1]*1.2)

fontsize('scale', 1.5)
end
