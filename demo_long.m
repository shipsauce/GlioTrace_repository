% GlioTrace longer demo for running ROI selection, stabilization, 
% cell tracking, classification, vessel segmentation, and 
% subsequent slice, cell and drug statistics + plot and visualize.
%
% @authors: Madeleine Skeppås
% @date: 31/10 2025

% Specify path to TIFs
file_path = 'path';

% If smaller regions of interest should be selected from a larger image,
% run select_stabilize and indicate whether you want to select
% regions manually or have a text-file of coordinates...
manual_selection = 1; % 1 - manual selection, 0 - coordinate file
region_size = 500;
coordinate_file = 'coordinates.txt';
output_path_stacks = '.../stabilized_stacks/';
stackfile = select_stabilize(file_path, output_path_stacks,...
                        manual_selection, region_size);

% ... otherwise, perform stabilization separately
stackfile = stabilize_tifs(file_path, output_path_stacks);

% Save stackfile for future use
save('stackfile.txt','stackfile');

% Define output path for generated videos (optional, otherwise run gliotrace without specifying path)
output_path_videos = '.../output/';

% Define metadata for all experiments
path_to_metadata = '../info_about_exp.xlsx';

% Run GlioTrace on stacks
video_output = '/Desktop/my_videos/';
detection_sensitivity = 0; % Tweak cell detection sensitivity in the tumor-channel (1 - most sensitive, 0 - least sensitive)
perturbations = "all";
gliotrace_output = gliotrace(stackfile', path_to_metadata, detection_sensitivity, perturbations, video_output);

% Plot example results
% Plot the relative proportions of cells in different classes for each celline
perturbation = {"control"};
style = "mean";
plot_proportions(gliotrace_output.slice_statistics,style,perturbation)

% Plot the speed + TAD of cells in different morphological classes for one celline
celline = "pat 1";
perturbation = "control";
dose = 0;
plot_cell_statistics(gliotrace_output.cell_statistics, celline, perturbation, dose)

% Plot drug effects
perturbation = "dasatinib";
doses = [16 32];
plot_drug_effects(gliotrace_output.drug_statistics, gliotrace_output.slice_statistics, perturbation, doses)


