function imout = customReadFunction_semseg(path)
    image = imread(path);
    imout = imresize(image, [250 250]);
end