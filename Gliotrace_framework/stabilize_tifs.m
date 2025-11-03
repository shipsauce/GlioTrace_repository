function stacktable = stabilize_tifs(file_path, output_path)
% @authors: Madeleine Skeppås
% @date: 31/10 2025

% Check whether folder contains TIF-files only or subfolders with TIFs
files = dir(fullfile(file_path, '**', '*.tif')); 
inSubfolder = any(~strcmp({files.folder}, file_path));

if(inSubfolder) %this case works
    folders = unique({files.folder});
    loop_length = length(folders);
else
    loop_length = 1;
    folders = file_path;
end

% Loop through TIF folder(s)
for k=1:loop_length
    files_sub = files(strcmp({files.folder}, folders{k}));
    foldername = regexp(files_sub(1).folder, '[^/\\]+$', 'match');
    experiment_of_interest = strrep(foldername{:}, ' ', '_');
 
    Tstack=[]; % tumor image stack
    Vstack=[]; % vasculature image stack
    Bstack=[]; % blue channel
    for time=1:length(files_sub)
        im_path=[files_sub(time).folder '/' files_sub(time).name];
        im=imread(im_path);
        Tstack(:,:,time)=im(:,:,2);
        Vstack(:,:,time)=im(:,:,1);
        Bstack(:,:,time)=im(:,:,3);
        time
    end

    % stabilize stack
    [Tstack_stab,Vstack_stab,Bstack_stab,~]=register_stack_DT(Tstack,Vstack,Bstack);
    
    stack = {};
    stack.Tstack = Tstack_stab;
    stack.Vstack = Vstack_stab;
    stack.Bstack = Bstack_stab;

    stackname = [output_path '/' experiment_of_interest '_roi_' num2str(k) '_stack.mat'];
    save(stackname, 'stack','-v7.3');
    stacktable(k) = cellstr(stackname);
end
end