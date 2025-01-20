function [vasc_length_stack, segstack] = quantify_vasculature(Vstack, sensitivity)
% This function takes the red channel of a brain slice stack and 
%
% Input parameters:
%   Vstack - red channel of stack
%   sensitivity - threshold for binarisation
%
% Output parameters:
%   vasc_length_stack - array of length measures of vasculature across
%   stack
%   segstack - the input stack with a blue overlay of the segmented vessels
%   
% @authors: Madeleine Skeppås
% @date: 17012025
    
% Toggle to enter debug mode
debug = 0;

vasc_length_stack = {};
segstack = [];

if(debug)
    tiledlayout(2,2,"TileSpacing","compact","Padding","compact")
end

% Iterate through the images of the stack
for i=1:size(Vstack,3)

    % Reduce noise and binarize with adaptive thresholding
    im = Vstack(:,:,i)/100;
    im = medfilt2(im);
    im = imgaussfilt(im,2);
    im_bin = imbinarize(im,"adaptive","ForegroundPolarity","bright","Sensitivity",sensitivity);

    % Remove objects that are too small or have to low eccentricity
    imlab = bwlabel(im_bin);
    regs = regionprops(imlab,'Eccentricity', 'Area');
    im_cleared = ismember(imlab, find((([regs.Eccentricity] > 0.5) .* [regs.Area] > 700)));
    
    % Skeletonize using the medial axis transform
    skel = bwskel(im_cleared,'MinBranchLength',10);
    
    % The length of the vasculature
    vasc_length_stack{i} = sum(sum(skel));

    % Create a grayscale image and the blue overlay
    imgray = im2gray(im);
    imgray = 0.5 * (imgray - 0.5) + 0.5;
    rgbImage = zeros([size(im) 3]);
    rgbImage(:,:,1) = imgray .* ~imdilate(skel,strel('disk',2)) + (65/255) * imdilate(skel,strel('disk',2));
    rgbImage(:,:,2) = imgray .* ~imdilate(skel,strel('disk',2)) + (105/255) * imdilate(skel,strel('disk',2));
    rgbImage(:,:,3) = imgray .* ~imdilate(skel,strel('disk',2)) + (225/255) * imdilate(skel,strel('disk',2));

    segstack(:,:,:,i) = im2uint8(rgbImage);
    
    % In debug mode the different steps are visualized
    if(debug)
        nexttile(1)
        imshowpair(im,skel)
        nexttile(2)
        imshowpair(im,im_cleared)
        nexttile(3)
        imshow(im)
        nexttile(4)
        plot(cell2mat(vasc_length_stack))
        hold on
    end
end

end