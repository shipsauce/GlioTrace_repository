function [count, imf, filtered_peaks, bw2] = count_cells(im, mask)
% This function takes an image and counts the number of fluorescent peaks.
%
% Input parameters:
% im - image wherein to count cells
% mask - binary mask for excluding/ including certain regions; specifically
%           for counting cells that express more than one marker
%
% Output parameters:
% count - number of cells found in image
% imf - filtered image
% filtered_peaks - image of identified peaks
% bw2 - binarized version of original image
%
% @authors: Madeleine Skeppås
% @date: 17012025

% Specify LoG filter corresponding to cell size
sigmah=3; 
hsizeh = sigmah*10; 
h = fspecial('log', hsizeh, sigmah);

% Gaussian filtering
imf=imgaussfilt(im,2);

% Mexican hat filter (LoG) for edge detection  
blob_im = conv2(imf,-h,'same');
bw = blob_im > imdilate(blob_im, [1 1 1; 1 0 1; 1 1 1]);

bw2 = imbinarize(imf,0.1);

% Apply mask if provided
if(nargin > 1)
    filtered_peaks = bw .* bw2 .* mask;
else
    filtered_peaks = bw .* bw2;
end

% Count number of detected peaks
count = sum(sum(filtered_peaks));

end