function vis_tracking_HD(traX,traY, mystack, myvasc, info, phenotypes, mode, start, path)
% This function takes the stack of time-lapse images and the associated cell
% coordinates and properties and creates a video visualizing the track of
% individual cells over time along with their morphological or TME
% association status at every time step.
% 
% Input parameters:
% traX - matrix with time-resolved cell x-coordinates (row=frame, col=cell)
% traY - matrix with time-resolved cell y-coordinates (row=frame, col=cell)
% mystack - green channel of stack
% myvasc - red channel of stack
% info - information on set, exp and ROI of stack; for writing to folder
% phenotypes - cell array of matrices of the time-resolved 
%                   properties of cells (row=frame, col=cell)
% mode - whether to visualize the morphological or TME classification
% start - starting index of stack; the first frame where cells appear
% path - output folder for videos
%
% @authors: Madeleine Skeppås
% @date: 14012025    

close all force
warning('off', 'MATLAB:audiovideo:VideoWriter:mp4FramePadded');

if(mode == "morphology")

    % Create a VideoWriter object and give it a filename for writing to
    % folder
    filename = [path '/set_' num2str(info.set) '_exp_' num2str(info.exp) '_roi_' num2str(info.roi) '_morphology.mp4'];
    obj=VideoWriter(filename,'MPEG-4');
    obj.FrameRate = 5;
    obj.Quality = 100;
    open(obj);
      
    % Color-coding classes (branching = red, diffuse translocation =
    % purple, junk = black, locomotion = blue, perivascular
    % translocation = green, round = yellow)

     colors = [
     1     0     0
     1     0     1
     0     0     0
     0     0     1
     0     1     0
     1     1     0];

    
    % Iterate through the stack
    for i=2:size(mystack,3)
        
        % Create an RGB image
        im = zeros([size(mystack,1,2) 3]);
        im(:,:,2)=medfilt2(mystack(:,:,i));
        im(:,:,1)=myvasc(:,:,i);
        im = uint8(im);

        % Check that we have passed the starting index
        if(i>=start)

            % Find index by correcting for the potential difference between
            % frame index and starting index of tracks
            shift = start-1;
            trax_i = i-shift;

            % Index the historic coordinates of cells (up to 5 frames back)
            ix=max(trax_i-5,1);
            ix=ix:trax_i;

            % Ignore tracks shorter than 4 frames
            idx = sum(~isnan(traX),1) > 3;
            traX_t = traX(:,idx);
            traY_t = traY(:,idx);
            labs = phenotypes{7}(:,idx);

            % Plot cell trajectories
            x=traX_t(ix,:);
            y=traY_t(ix,:);

            % Draw the blue lines on the image
            for k = 1:size(x, 2) - 1
                coords = [x(~isnan(x(:,k)),k) y(~isnan(y(:,k)),k)];
                if(size(coords,1) > 1)
                    % Add the blue lines to the snippet using insertShape
                    im = insertShape(im, 'Line', [x(~isnan(x(:,k)),k) y(~isnan(y(:,k)),k)], 'Color', 'blue', 'LineWidth', 2);
                end
            end
            
            % Find the associated morphological labels of cells for the
            % current frame
            ind = labs(trax_i,:);
            clrs = colors(ind(~isnan(ind)),:);
            xs = traX_t(trax_i,:);
            ys = traY_t(trax_i,:);

            x = xs(~isnan(xs));
            y = ys(~isnan(ys));
            
            % Plot the cell morphological labels as a scatter of different
            % colored spots
            im = insertShape(im, 'filled-circle', [x' y' repmat(3, length(x),1)], 'Color', clrs*255, 'Opacity', 1);
        end
        
        curr_frame = im2frame(imresize(im, 'scale', 10));
        writeVideo(obj,curr_frame);
    end

    close(obj);

elseif(mode == "tme")
    % Create a VideoWriter object and give it a filename for writing to
    % folder
    filename = [path '/set_' num2str(info.set) '_exp_' num2str(info.exp) '_roi_' num2str(info.roi) '_tme.mp4'];
    obj=VideoWriter(filename,'MPEG-4');
    obj.FrameRate = 5;
    obj.Quality = 100;
    open(obj);

    % Color-coding classes (microglia colocalized = yellow, vessel
    % associated = green, non-associated = white)
    colors = [
    1     1     0
    0     1     0
    1     1     1];  
    
    % Iterate through the stack
    for i=2:size(mystack,3)

        % Plot TME labels on an RGB image
        im2 = zeros([size(mystack,1,2) 3]);
        im2(:,:,2)=medfilt2(mystack(:,:,i));
        im2(:,:,1)=myvasc(:,:,i);
        im2 = uint8(im2);

        if(i>=start)
            shift = start-1;
            trax_i = i-shift;

            % Ignore tracks shorter than 4 frames
            idx = sum(~isnan(traX),1) > 3;
            traX_t = traX(:,idx);
            traY_t = traY(:,idx);
            labs = phenotypes{8}(:,idx);
            
            ind = labs(trax_i,:);

            clrs = colors(ind(~isnan(ind)),:);

            xs = traX_t(trax_i,:);
            ys = traY_t(trax_i,:);

            x = xs(~isnan(xs));
            y = ys(~isnan(ys));
            
            % Plot the cell morphological labels as a scatter of different
            % colored spots
            im2 = insertShape(im2, 'circle', [x' y' repmat(7, length(x),1)], 'Color', clrs*255, 'Opacity', 1, 'LineWidth', 3);
        end

        % Plot TME labels on the red channel only
        im3 = zeros([size(mystack,1,2) 3]);
        im3(:,:,1)=myvasc(:,:,i);
        im3 = uint8(im3);

        if(i>=start)
            shift = start-1;
            trax_i = i-shift;
    
             % Ignore tracks shorter than 4 frames
            idx = sum(~isnan(traX),1) > 3;
            traX_t = traX(:,idx);
            traY_t = traY(:,idx);
            labs = phenotypes{8}(:,idx);
            
            ind = labs(trax_i,:);

            clrs = colors(ind(~isnan(ind)),:);

            xs = traX_t(trax_i,:);
            ys = traY_t(trax_i,:);

            x = xs(~isnan(xs));
            y = ys(~isnan(ys));
            
            % Plot the cell morphological labels as a scatter of different
            % colored spots
            im3 = insertShape(im3, 'circle', [x' y' repmat(7, length(x),1)], 'Color', clrs*255, 'Opacity', 1, 'LineWidth', 3);
        end
        
        im_merged=imtile({im2, im3},'bordersize',[3 3],'GridSize',[1,2],'backgroundcolor','white');
        curr_frame = im2frame(im_merged);
        writeVideo(obj,curr_frame);
    end
    close(obj);
end
end