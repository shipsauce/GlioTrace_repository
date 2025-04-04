function [original_sequence_paths, stable_sequence_paths, ROI_sequence_paths, segmented_ROI_sequence_paths]=get_paths_to_data(metadata)
% This function takes a path to an excel file containing metadata as input,
% and returns the paths of the images.
%
% @authors: Sven Nelander
% @date: 140624

original_sequence_paths=metadata.path;%textread(file_with_exp_paths,'%s\n','delimiter','\t');

base='/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/';
n=length(original_sequence_paths);
for i=1:n
    original_sequence_paths{i}=[base original_sequence_paths{i} '/'];
end

stable_sequence_paths={};
for i=1:n
    stable_sequence_paths{i,1}=[ '/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/stabilized/exp' num2str(i) '/'];
end

ROI_sequence_paths={};
for i=1:n
    ROI_sequence_paths{i}=[ '/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/ROI/exp' num2str(i) '/'];
end

segmented_ROI_sequence_paths={};
for i=1:n
    segmented_ROI_sequence_paths{i}=[ '/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/segmented/exp' num2str(i) '/'];
end
