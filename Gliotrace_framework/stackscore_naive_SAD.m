function sad=stackscore_naive_SAD(gbm)
% This function calculates the mean sum of absolute differences (SAD)
% across the full stack.
%
% Input parameters:
% gbm - image stack
%
% Output parameters:
% sad - scalar value of mean SAD
%
% @authors: Sven Nelander
% @date: 14082024

    % Filter away low intensity pixels and set to zero
    gbm=gbm-20;
    gbm(gbm<0)=0;
    
    % Create an image stack of absolute differences between each pair of
    % images i and i+1
    sad=abs(gbm(:,:,2:end)-gbm(:,:,1:end-1));
    
    % Calculate the mean SAD across the entire stack
    sad=mean(mean(mean(sad,1),2),3);
end

