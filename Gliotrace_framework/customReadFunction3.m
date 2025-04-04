function imout = customReadFunction3(path)
    load("mean_image_tme.mat")
    image = imread(path);
    sz = size(image);

    if(sz(1) < 61)
        imout = im2double(imresize(image, [61 61]));
        imout = imout - mean_image;
    else
        imout = im2double(image);
        imout = imout - mean_image;
    end
end