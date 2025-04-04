function stabilize_IMX_data(mypath,outpath, crop_rectangle,working_scale)
% This function stbilizes images of a stack and writes the resulting video
% as a tif into a specified output folder.
%
% @authors: Sven Nelander
% @date: 140624

files_in=getfilenames(mypath);


image=imread([mypath files_in{1}]);
if(length(crop_rectangle))
    image=imcrop(image,crop_rectangle);
end
image=imresize(image,'Scale',working_scale);
imwrite(image,[outpath files_in{1}],'tiff');

IMred=image(:,:,1);
% Process all frames in the video
movMean = double(IMred(:,:,1))/255*2;
imgB = movMean;
imgBp = imgB;
correctedMean = imgBp;
ii = 2;
Hcumulative = eye(3);
IMredstable=[];
IMgreenstable=[];
while ii<length(files_in)
    % Read in new frame
    imgA = imgB; % z^-1
    imgAp = imgBp; % z^-1
    
    try
        image=imread([mypath files_in{ii}]);
        
        if(length(crop_rectangle))
            image=imcrop(image,crop_rectangle);
        end
        image=imresize(image,'Scale',working_scale);
        IMred=image(:,:,1);
        IMgreen=image(:,:,2);
        
        imgB = double(IMred)/255*2;
        movMean = movMean + imgB;
        
        % Estimate transform from frame A to frame B, and fit as an s-R-t
        if(ii==31)
            1
        end
        [H,HsRt] = cvexEstStabilizationTform_modified(imgA,imgB,0.1);
        %HsRt = cvexTformToSRT(H);
        Hcumulative = HsRt * Hcumulative;
        imgBp = imwarp(imgB,affine2d(Hcumulative),'OutputView',imref2d(size(imgB)));
        Hcumulative
        % Display as color composite with last corrected frame
        imshow(imfuse(imgAp,imgBp,'ColorChannels','red-cyan'));
        IMredstable=imgBp;
        IMgreenstable=imwarp(double(IMgreen)/255*2,affine2d(Hcumulative),'OutputView',imref2d(size(imgB)));
        
        IM=IMredstable;
        IM(:,:,2)=IMgreenstable;
        IM(:,:,3)=zeros;
        IM=uint8(IM*255/2);
        imwrite(IM,[outpath files_in{ii}],'tiff');
        
        drawnow;
        correctedMean = correctedMean + imgBp;
    catch
        warning('something went wrong');
    end
    'frame number'
    ii = ii+1
end
correctedMean = correctedMean/(ii-2);
movMean = movMean/(ii-2);

end