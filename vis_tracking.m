function vis_tracking(traX,traY, mystack, myvasc, info, phenotypes, mode, start)
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
%
% @authors: Madeleine Skeppås
% @date: 14012025    

close all force

if(mode == "morphology")

    % Create a VideoWriter object and give it a filename for writing to
    % folder
    folder = ['/Users/madsk418/UU Dropbox/Madeleine S/Simulation_and_invasion/comp/output/Madeleine/Tracking_vids/Single_cell_revision/morph/'];
    filename = [folder 'set_' num2str(info.set) '_exp_' num2str(info.exp) '_roi_' num2str(info.roi) '.mp4'];
    obj=VideoWriter(filename,'MPEG-4');
    obj.FrameRate = 5;
    open(obj);
    figure('Renderer', 'painters', 'Position', [500 500 1000 800])
    sgtitle(['Set ' num2str(info.set) ' exp ' num2str(info.exp) ' roi ' num2str(info.roi)])
     
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
        im(:,:,2)=medfilt2(mystack(:,:,i))/100;
        im(:,:,1)=myvasc(:,:,i)/100;
        imshow(im);
        hold on;

        % Check that we have passed the starting index
        if(i>=start)

            % Find index by correcting for the potential difference between
            % frame index and starting index of tracks
            shift = start-1;
            trax_i = i-shift;

            % Index the historic coordinates of cells (up to 5 frames back)
            ix=max(trax_i-5,1);
            ix=ix:trax_i;

            % Plot cell trajectories
            plot(traX(ix,:),traY(ix,:),'b-','LineWidth',3);
            
            % Find the associated morphological labels of cells for the
            % current frame
            ind = phenotypes{7}(trax_i,:);
            clrs = colors(ind(~isnan(ind)),:);
            xs = traX(trax_i,:);
            ys = traY(trax_i,:);

            % Plot the cell morphological labels as a scatter of different
            % colored spots
            scatter(xs(~isnan(xs)),ys(~isnan(ys)),50,clrs,'filled', 'MarkerEdgeColor','k')
        end
        
        drawnow;
        frame=getframe(gcf);
        writeVideo(obj,frame);
    end

    close(obj);

elseif(mode == "tme")
    folder = ['/Users/madsk418/UU Dropbox/Madeleine S/Simulation_and_invasion/comp/output/Madeleine/Tracking_vids/Final_TME_inspection_jan_2025/' info.HGCC{1} '/' info.perturbation{:} '/'];
    filename = [folder 'set_' num2str(info.set) '_exp_' num2str(info.exp) '_roi_' num2str(info.roi) '.mp4']; 
    obj=VideoWriter(filename,'MPEG-4');
    obj.FrameRate = 5;
    open(obj);
    figure('Renderer', 'painters', 'Position', [500 500 1000 800])

    % Color-coding classes (microglia colocalized = yellow, vessel
    % associated = green, non-associated = white)
    colors = [
    1     1     0
    0     1     0
    1     1     1];  
    
    tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    sgtitle(['Set ' num2str(info.set) ' exp ' num2str(info.exp) ' roi ' num2str(info.roi)])
    
    % Iterate through the stack
    for i=2:size(mystack,3)

        % Plot TME labels on an RGB image
        nexttile(1);
        im = zeros([size(mystack,1,2) 3]);
        im(:,:,2)=medfilt2(mystack(:,:,i))/100;
        im(:,:,1)=myvasc(:,:,i)/100;
        imshow(im);
        hold on;

        if(i>=start)
            shift = start-1;
            trax_i = i-shift;
            
            ind = phenotypes{8}(trax_i,:);

            clrs = colors(ind(~isnan(ind)),:);

            xs = traX(trax_i,:);
            ys = traY(trax_i,:);
            scatter(xs(~isnan(xs)),ys(~isnan(ys)),500,clrs, LineWidth=3)
        end

        % Plot TME labels on the red channel only
        nexttile(2);
        im = zeros([size(mystack,1,2) 3]);
        im(:,:,1)=myvasc(:,:,i)/100;
        imshow(im);
        hold on;

        if(i>=start)
            shift = start-1;
            trax_i = i-shift;
    
            ind = phenotypes{8}(trax_i,:);
    
            clrs = colors(ind(~isnan(ind)),:);
    
            xs = traX(trax_i,:);
            ys = traY(trax_i,:);
            scatter(xs(~isnan(xs)),ys(~isnan(ys)),500,clrs, LineWidth=3)
        end
        
        drawnow;

        frame=getframe(gcf);
        writeVideo(obj,frame);
    end
    close(obj);
end
end
