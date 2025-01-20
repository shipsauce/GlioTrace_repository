function [cellsx,cellsy,intensity,FEAT,VASC]=macro_track2(mystack,vasc,sigmah,hsizeh,cutoff, blocksize, mode)
% This function takes the green channel of a brain slice stack and
% identifies cell coordinates of cell bodies. 
% 
% Input parameters:
% mystack - green channel of stack
% vasc - red channel of stack
% sigmah - standard dev. of LoG filter
% hsizeh - size of LoG filter
% cutoff - intensity threshold for cell detection
% blocksize - the size of the extracted snapshot from extractFeatures
% mode - 'normal' or 'sparse'; when cells are unusually faint and/or sparse
% in the video, an alternative set of parameters are used for cell
% detection
%
% Output parameters:
% cellsx - cell array with frame-wise vectors of x-coordinates
% cellsy - cell array with frame-wise vectors of y-coordinates
% intensity - pixel-values at x,y-coordinates
% FEAT - cell array storing the extracted snapshots as one matrix per frame
%           where each row contains the compressed snapshot of one cell
% VASC - same as FEAT but snapshot is extracted from the red channel
%
% In debug mode, every frame is displayed as green channel only + green and
% red channel combined, and every detected cell is marked by a white circle
%
% @authors: Madeleine Skeppås, Sven Nelander
% @date: 10012025    

% Toggle to enter debug mode
debug=false;    

% Create a Logarithm-of-Gaussian filter
    h1 = fspecial('log', hsizeh, sigmah);

% Create structures for storing cell coordinates and associated intensities
cellsy={};
cellsx={};
intensity={};

% Iterate through the frames of the stack
for i=1:size(mystack,3)

    % Create temporary variables for storing cell coordinates
    xs = [];
    ys = [];
    
    im=mystack(:,:,i); % Green channel
    im_vasc=vasc(:,:,i); % Red channel

    % Apply LoG filter
    blob_im_1 = conv2(double(im)/255,-h1,'same');
  
    % Blob detection operation
    bw = blob_im_1 > imdilate(blob_im_1, [1 1 1; 1 0 1; 1 1 1]);

    % Create binary masks identifying cell bodies and ignoring protrusions
    cellbodypixels_1=imdilate(imerode(blob_im_1>cutoff,strel('disk',3)),strel("disk",3));
    cellbodypixels_2 = imbinarize(im/255);
    cellbodypixels_3 = imdilate(cellbodypixels_2, strel('disk',3));
    
    % If cells are sparse in frames, use another mask for cell bodies
    if(mode == "sparse")
        cellbodypixels = cellbodypixels_1;
    else
        cellbodypixels = cellbodypixels_3;
    end

    % Label cell bodies as individual objects
    L=bwlabel(cellbodypixels);

    % Measure properties of individual objects
    regs = regionprops('table',L, 'Area');

    % Ignore potential cell bodies with an area above some threshold
    idx = find(regs.Area > 3000);
    if(~isempty(idx))
        for j=1:length(idx)
            L(L==idx(j)) = 0;
        end
    end

    cellbodypixels = L;

    % Perform Gaussian filtering of green channel to identify bulk regions
    imsmooth = imgaussfilt(im,30);
    mask = imbinarize(uint8(imsmooth)); % Mask to ignore bulk regions

    % If max intensity of Gaussian filtered image exceeds threshold, this
    % is indicative of a bulk region, and the bulk mask should be applied
    % to prune away coordinates in this region
    if(max(max(imsmooth)) > 70)
        % Iterate through the cell body objects
        for k=1:max(max(cellbodypixels))
            px = cellbodypixels == k;

            % Prune coordinates using the cell body, the binarized original
            % image (cellbodypixels_2) and the bulk mask
            [y,x] = find((bw).*(blob_im_1>cutoff).*px.*cellbodypixels_2.*imcomplement(mask));
            
            % If still more than one coordinate for the cell body object
            if(length(x)>1)
                xmean=round(mean(x));
                ymean=round(mean(y));

                area2 = sum(sum(px));
                
                % If the mean coordinate is present within the cell body
                % and the area of the object is within some thresholds,
                % keep the mean coordinate instead
                if(cellbodypixels_3(ymean,xmean) && area2 > 200 && area2 < 1000)
                    xs = [xs; xmean];
                    ys = [ys; ymean];
                else % If the mean coordinate is outside of the cell body, 
                    % or the area of the object is outside of area boundaries,
                    % this indicates that there could be more than one cell
                    % hiding within the cell body object. In that case,
                    % keep the coordinates as is.
                    xs = [xs; x];
                    ys = [ys; y];
                end
            else
                try
                xs = [xs; x];
                ys = [ys; y];
                catch
                    1;
                end
            end
        end

    else
        % Iterate through the cell body objects
         for k=1:max(max(cellbodypixels))
            px = cellbodypixels == k;

            % Prune coordinates using the cell body, the binarized original
            % image (cellbodypixels_2) and the bulk mask
            [y,x] = find((bw).*(blob_im_1>cutoff).*px.*cellbodypixels_2);

            % If still more than one coordinate for the cell body object
            if(length(x)>1)
                xmean=round(mean(x));
                ymean=round(mean(y));

                area2 = sum(sum(px));
                
                if(cellbodypixels_3(ymean,xmean) && area2 > 200 && area2 < 1000)
                    xs = [xs; xmean];
                    ys = [ys; ymean];
                else
                    xs = [xs; x];
                    ys = [ys; y];
                end

            else
                try
                xs = [xs; x];
                ys = [ys; y];
                catch
                    1;
                end
            end
         end
    end
    
    % As long as frame is not empty on cells, extract snapshots
    if(~isempty(cellsy{i}))
        [featuresA, pointsA] = extractFeatures(medfilt2(im), [xs ys],'Method','Block','BlockSize',blocksize);
        [featuresB, pointsB] = extractFeatures(medfilt2(im_vasc), [xs ys],'Method','Block','BlockSize',blocksize);
        
        % Save coordinates as points from extractFeatures in case some
        % coordinates have been excluded from snapshot extraction because
        % they were too close to the edge
        cellsx{i}=pointsA(:,1);
        cellsy{i}=pointsA(:,2);
        
        % Save snapshots and intensity values of coordinates
        FEAT{i}=featuresA;
        VASC{i}=featuresB;
        f2 = sub2ind(size(im2),pointsA(:,1), pointsB(:,2));
        intensity{i} = im2(f2);
    end
    
    % Debug mode displays green channel separately & green and red channel
    % combined
    if(debug)
        subplot(1,2,1)
        jm=cat(3,medfilt2(im_vasc),medfilt2(im),0*im);
        imshow(jm/255);
        hold on;
        scatter(cellsx{i}, cellsy{i}, 'wo')
        hold off;
        subplot(1,2,2);
        green = zeros([size(im) 3]);
        green(:,:,2) = medfilt2(im);
        imshow(green/255)
        hold on;
        scatter(cellsx{i}, cellsy{i}, 'wo')
        hold off;
        drawnow;
        pause(0.4)
    end
end

% If no cells could be detected in the stack, throw this error and restart
% the algorithm with another set of parameters
if(all(cellfun(@(x) isequal(x, zeros(0, 1)) || isequal(x, zeros(0,0)) , cellsx)))
    error('Output argument "cellsx" (and possibly others) not assigned a value in the execution with "macro_track2" function.')
end
end