function [vasc_length_stack, segstack] = quantify_vasculature_v3(Vstack, info, net_trained)
% This function takes the red channel of a brain slice stack and 
%
% Input parameters:
%   Vstack - red channel of stack
%   info - threshold for binarisation
%   net_trained - 
%
% Output parameters:
%   vasc_length_stack - array of length measures of vasculature across
%   stack
%   segstack - the input stack with a blue overlay of the segmented vessels
%   
% @authors: Madeleine Skeppås
% @date: 15022025

debug = 1;

vasc_length_stack = {};
segstack = [];
close all force

if(debug)
    close all force
    figure('Renderer', 'painters', 'Position', [500 500 1500 1500])
    tiledlayout(2,2,"TileSpacing","compact","Padding","compact")
    folder = ['/Users/madsk418/UU Dropbox/Madeleine S/Simulation_and_invasion/comp/output/Madeleine/Tracking_vids/Final_vasc_inspection_feb_2025/' info.HGCC{1} '/'];
    filename = [folder 'set_' num2str(info.set) '_exp_' num2str(info.exp) '_roi_' num2str(info.roi) '_' info.vessel_origin{:} '.mp4'];
    obj=VideoWriter(filename,'MPEG-4');
    obj.FrameRate = 5;
    open(obj);
end

% Iterate through the images of the stack
for i=1:size(Vstack,3)

    % Reduce noise and binarize with adaptive thresholding
    im = medfilt2(Vstack(:,:,i));

    [h, w, c] = size(im);
    halfH = floor(h/2);
    halfW = floor(w/2);

    imgParts = {im(1:halfH, 1:halfW, :), im(1:halfH, halfW+1:end, :), ...
                im(halfH+1:end, 1:halfW, :), im(halfH+1:end, halfW+1:end, :)};
    
    predictedMask = semanticseg(uint8(im), net_trained);
    predictedMask = ~(predictedMask == 'C1');
    predictedMask = imclose(predictedMask, strel('disk',3));
    
    % Skeletonize using the medial axis transform
    skel = bwskel(predictedMask,'MinBranchLength', 15);
    
    % The length of the vasculature
    vasc_length_stack{i} = sum(sum(skel));

    % Create a grayscale image and the blue overlay
    imgray = im2gray(im/100);
    imgray = 0.5 * (imgray - 0.5) + 0.5;
    skel_resized = imresize(skel, [size(imgray)]);
    rgbImage = zeros([size(im) 3]);
    rgbImage(:,:,1) = imgray .* ~imdilate(skel_resized,strel('disk',2)) + (65/255) * imdilate(skel_resized,strel('disk',2));
    rgbImage(:,:,2) = imgray .* ~imdilate(skel_resized,strel('disk',2)) + (105/255) * imdilate(skel_resized,strel('disk',2));
    rgbImage(:,:,3) = imgray .* ~imdilate(skel_resized,strel('disk',2)) + (225/255) * imdilate(skel_resized,strel('disk',2));

    segstack(:,:,:,i) = im2uint8(rgbImage);
    
    % In debug mode the different steps are visualized
    if(debug)
        nexttile(1)
        imshowpair(im,skel)
        nexttile(2)
        imshowpair(im,predictedMask)
        nexttile(3)
        imshow(cat(3, im, zeros(size(im)), zeros(size(im)))./255)
        nexttile(4)
        if(i >= 5)
            baseline = median(cell2mat(vasc_length_stack(1:5)));
            plot((cell2mat(vasc_length_stack) ./ baseline) .* 100)
        end
        ylim([0 150])
        hold on
        drawnow
        frame=getframe(gcf);
        writeVideo(obj,frame);
    end
    i
end

if(debug)
    close(obj);
end

end