function [vasc_length_stack, segstack] = segment_quantify_vasculature(Vstack, info, path)
% This function takes the red channel of a brain slice stack and uses a
% neural network to segment the vasculature with semantic segmentation. The
% resulting binary mask is measured in length by skeletonization followed
% by summation of non-zero pixels. The binary mask is also used to create
% an overlay to the red channel, which is saved as a a video into the
% specified path. Both the vasculature length across the time lapse and the
% overlay image stack is saved into a table which is specified as output.
%
% Input parameters:
%   Vstack - red channel of stack
%   info - metadata for the video title
%   path - output path of generated video
%
% Output parameters:
%   vasc_tbl with variables:
%       vasc_length_stack - array of length measures of vasculature across
%                           stack
%       segmented_stack - the input stack with a blue overlay of the segmented vessels
%   
% @authors: Madeleine Skeppås
% @date: 15022025

% Load the neural network for semantic segmentation of vessels
load("net_2025_02_20.mat")

vasc_length_stack = {};
segstack = [];
close all force

if(~isempty(path))
    filename = [path '/set_' num2str(info.set) '_exp_' num2str(info.exp) '_roi_' num2str(info.roi) '_vasculature_segmentation.mp4'];
    obj=VideoWriter(filename,'MPEG-4');
    obj.FrameRate = 5;
    obj.Quality = 50;
    open(obj);
end

% Iterate through the images of the stack
for i=1:size(Vstack,3)

    im = medfilt2(Vstack(:,:,i)); % Median filtering

    % Predict vasculature mask using segmantic segmentation with a U-net
    predictedMask = semanticseg(uint8(im), net_trained);
    predictedMask = ~(predictedMask == 'C1');

    % Post-process the mask with morphological closing
    predictedMask = imclose(predictedMask, strel('disk',3));
    
    % Skeletonize using the medial axis transform
    skel = bwskel(predictedMask,'MinBranchLength', 15);
    
    % The length of the vasculature is the sum of all non-zero pixels
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
    
    if(~isempty(path))
        % Write to video
        im1 = zeros([size(im) 3]);
        im1(:,:,1) = im;
        im1 = uint8(im1);
    
        im_merged=imtile({im1, rgbImage},'bordersize',[3 3],'GridSize',[1,2],'backgroundcolor','white');
        curr_frame = im2frame(im_merged);
        writeVideo(obj,curr_frame);
    end

end

if(~isempty(path))
    close(obj);
end


end