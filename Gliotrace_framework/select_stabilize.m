function stacktable = select_stabilize(file_path, output_path, manual_selection, region_size, coordinate_file)
% Region selection and stabilization. 
%
% Function can handle TIFs organized in two ways;
%   1 - file_path is one folder with TIFs from the same experiment
%   2 - file_path contains subfolders with TIFs from separate experiments
%
% Folder name becomes experiment_id in stacktable.
% 
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

% Handle last input argument
if(~manual_selection)
    coordinate_file = readtable(region_size);
end

stacktable_total = [];
ROI_table = table();

% Loop through TIF folder(s)
for i=1:loop_length
    
    % Manual region selection
    if(manual_selection)
        files_sub = files(strcmp({files.folder}, folders{i}));
        time=min(50,length(files_sub));
        im=imread([files_sub(time).folder '/' files_sub(time).name]);
        foldername = regexp(files_sub(time).folder, '[^/\\]+$', 'match');
        experiment_of_interest = strrep(foldername{:}, ' ', '_');
        myROI = table();
        update=1;
        roi_count=1;
        stacktable = {};
       
        while(update)
            imshow(2*im)
            set(gcf, 'Windowstyle', 'normal')
            h = images.roi.Rectangle(gca,'Position',[50 50 region_size region_size],'StripeColor','r');
            update = 0;
    
            fprintf('Select one of the below and press enter\n');
            r=input('1=one more, 2=review all, 9=done\n');
        
            if(r==1)
                myROI.experiment{roi_count} = cellstr(experiment_of_interest);
                myROI{roi_count, {'X','Y','W','H'}} = num2cell(h.Position);
                roi_count = roi_count+1;
                update=1;
            end
        
            if(r==2)
                myROI.experiment{roi_count} = cellstr(experiment_of_interest);
                myROI{roi_count, {'X','Y','W','H'}} = num2cell(h.Position);
                roi_count = roi_count+1;
                update=1;

                imshow(2*im)
                set(gcf, 'Windowstyle', 'normal')

                for m=1:height(myROI)
                    h(m) = images.roi.Rectangle(gca,'Position', cell2mat(myROI{m,2:end}),'StripeColor','r');
                end

                fprintf('Regions adjusted y/n?\n');
                r=input('1=yes, 2=no\n');
                
                if(r==1)
                    for m=1:height(myROI)
                        myROI{roi_count, {'X','Y','W','H'}} = num2cell(h(m).Position);
                    end
                end
            end
        
            if(r==9)
                myROI.experiment(roi_count) = cellstr(experiment_of_interest);
                myROI{roi_count, {'X','Y','W','H'}} = num2cell(h.Position);
                roi_count = roi_count+1;
                update=0;
            end
        end
    
        ROI_table=[ROI_table; myROI];

        for k=1:size(myROI,1)
            ROI = myROI(k,2:end);
            
            Tstack=[]; % tumor image stack
            Vstack=[]; % vasculature image stack
            Bstack=[]; % blue channel
            for time=1:length(files_sub)
                im_path=[files_sub(time).folder '/' files_sub(time).name];
                im=imread(im_path,'PixelRegion',{[ROI.Y ROI.Y+ROI.H],[ROI.X ROI.X+ROI.W]});
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

        stacktable_total = [stacktable_total stacktable];

    % Automatic region selection using coordinate file    
    else
        files_sub = files(strcmp({files.folder}, folders{i}));
        foldername = regexp(files_sub(1).folder, '[^/\\]+$', 'match');
        experiment_of_interest = strrep(foldername{:}, ' ', '_');

        ROIs_sub = coordinate_file(strcmp(coordinate_file.exp, experiment_of_interest),:);

        for k=1:height(ROIs_sub)
            ROI = ROIs_sub(k,2:end);
            
            Tstack=[]; % tumor image stack
            Vstack=[]; % vasculature image stack
            Bstack=[]; % blue channel
            for time=1:length(files_sub)
                im_path=[files_sub(time).folder '/' files_sub(time).name];
                im=imread(im_path,'PixelRegion',{[ROI.Y ROI.Y+ROI.H],[ROI.X ROI.X+ROI.W]});
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

        stacktable_total = [stacktable_total stacktable];
    end
end

end