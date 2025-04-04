function b=stackscore_naive_growth_rate(gbm,dt)
% This function calculates the growth rate across the image stack based on
% the mean intensity in each image.
%
% Input parameters:
% gbm - image stack
%
% Output parameters:
% sad - scalar value of mean SAD
%
% @authors: Sven Nelander
% @date: 14082024
    
    y=(squeeze(mean(mean(gbm,1),2)));
    t=[0:dt:dt*(size(gbm,3)-1)]';
    beta=pinv([ones(size(t)) t])*y;
    b=beta(2);
end
