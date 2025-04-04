function visualize_stack(gbm, vstack)
% Given one or two image stacks, this function visualizes the (combined)
% time-lapse series of images.
% 
% Input parameters: 
% - gbm - green channel
% - vstack - red channel
%
% @authors: Madeleine Skeppås
% @date: 20022025

for i=1:size(gbm,3)

    if(nargin==2)
        im = zeros([size(gbm,1,2) 3]);
        im(:,:,1) = vstack(:,:,i)/255;
        im(:,:,2) = gbm(:,:,i)/255;
    else
        im = zeros([size(gbm,1,2) 3]);
        im(:,:,1) = gbm(:,:,i)/255;
    end

    imshow(im)
    hold on
    drawnow
    pause(0.1)
end

end