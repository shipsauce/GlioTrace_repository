function imout = customReadFunction(path)
    image = imread(path);
    sz = size(image);

    if(sz(1) < 61)
        imout = im2double(imresize(image, [61 61]));
    else
        imout = im2double(image);
    end
end