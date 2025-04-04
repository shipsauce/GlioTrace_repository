function [Tstack_stab, Vstack_stab] = save_stack(roi, output_path, duplicates)
% This function takes a file of ROI coordinates and an output path to
% write the resulting matlab object to. From the full image, it crops
% smaller regions which is stabilized using register_stack and then written
% to the output path. If duplicates is set to true, it also checks whether
% any of the ROIs in the ROI-file have already been written to the output
% path.
%
% @authors: Sven Nelander, Madeleine Skeppås
% @date: 140624

if(nargin < 3)
duplicates = false;
end

metadata=readtable('/Volumes/MyGroups$/Iron/konfokalmikroskop/Hitesh Montage and Overlays/hitesh_metadata.xlsx');

fileList = dir(fullfile(output_path, '*.mat'));

[original_sequence_paths, stable_sequence_paths]=get_paths_to_data(metadata);

ROI_table = readtable(roi);

if(duplicates)
    fileList = dir(fullfile(output_path, '*.mat'));
    
    oldExps = [];
    
    for j=1:length(fileList)
        oldExps(j) = str2num(extractBetween(string(fileList(j).name), 'exp_', '_roi'));
    end
    
    exps = unique(ROI_table.exp);
    oldExps = unique(oldExps');
   
    newExps = exps(~ismember(exps,oldExps));

    ROI_table = ROI_table(ismember(ROI_table.exp,newExps),:);
end

for i=1:size(ROI_table,1)
    exp = ROI_table{i,1};
    files=getfilenames(stable_sequence_paths{exp});
    ROI = ROI_table(i,2:end);
    
    % load the stack (takes a minute or so)
    Tstack=[]; % tumor image stack
    Vstack=[]; % vasculature image stack
    for time=1:length(files)
        myimage=[stable_sequence_paths{exp} files{time}];
        im=imread(myimage,'PixelRegion',{[ROI.Y ROI.Y+ROI.H],[ROI.X ROI.X+ROI.W]});
        Tstack(:,:,time)=im(:,:,2);
        Vstack(:,:,time)=im(:,:,1);
        time
    end

    % stabilize stack
    [Tstack_stab,Vstack_stab,delta]=register_stack(Tstack,Vstack);
    
    stack = {};
    stack.Tstack = Tstack_stab;
    stack.Vstack = Vstack_stab;

    save([output_path '/' 'exp_' num2str(exp) '_roi_' num2str(i) '_stack.mat'], 'stack','-v7.3');

    i
end

end