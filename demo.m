% GlioTrace demo for running cell tracking, classification, vessel
% segmentation, and subsequent slice, cell and drug statistics.
%
% @authors: Madeleine Skeppås
% @date: 15022025

% Define paths to image stacks
stackfile = 'saved_stacks_complete.txt';

% Define output path for generated videos (optional, otherwise run gliotrace without specifying path)
output_path = '/Users/madsk418/Desktop/output/';

% Run GlioTrace on stacks
gliotrace_output = gliotrace(stackfile);

%%

% Plot the relative proportions of cells in different classes for each celline
perturbation = {"control","thapsigargin"};
style = "mean";
plot_proportions(gliotrace_output.slice_statistics,style,perturbation)

% Plot the speed + TAD of cells in different morphological classes for one celline
celline = "U3013MG";
perturbation = "control";
dose = 0;
plot_cell_statistics(gliotrace_output.cell_statistics, celline, perturbation, dose)

% Plot drug effects
perturbation = "dasatinib";
doses = [16 32];
plot_drug_effects(gliotrace_output.drug_statistics, gliotrace_output.slice_statistics, perturbation, doses)
