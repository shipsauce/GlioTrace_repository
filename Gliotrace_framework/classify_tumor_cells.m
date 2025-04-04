function properties = classify_tumor_cells(feat, vasc, blocksize, morph_net, trainedNetwork_tme, m)
% This function takes the extracted snapshots from macro_track2 and
% classifies them two times, once based on morphology and once based on
% the relation to the surrounding tumor microenvironment as represented by
% the stained vasculature in the red channel.
% 
% Input parameters:
% feat - cell array storing the extracted snapshots as one matrix per frame
%           where each row contains the compressed snapshot of one cell
% vasc - same as feat but snapshot is extracted from the red channel
% blocksize - size of snapshot
% morph_net - convolutional neural network for classifying morphology
% trainedNetwork_tme - convolutional neural network for classifying the
%                       relation to the tumor microenvironment
% m - stack index; used when saving snapshots for manual inspection
%
% Output parameters:
% properties - cell array where every cell contains a matrix with the
%               associated features of the cells in that frame, where rows
%               are cells and cols are features (class probabilities and
%               labels)
%
% In debug mode the classified snapshot and its associated morphology label
% is displayed for inspection. 

% @authors: Madeleine Skeppås
% @date: 13012025    

% Toggle to enter debug mode
debug = false;

% Load mean images for mean subtraction centering before classification
mean_image_6class = load('mean_image_6class.mat');
mean_image_6class = mean_image_6class.mean_image;
mean_image_tme = load('mean_image_tme.mat');
mean_image_tme = mean_image_tme.mean_image;

% Path used when saving snapshots
output_path = '/Users/madsk418/UU Dropbox/Madeleine S/Simulation_and_invasion/comp/output/Madeleine/Networks_validation_v5/';

properties={}; % Structure for saving classification results

% Iterate through the features of each frame
for i=1:length(feat)

    mat = feat{i}; % Snapshot from green channel
    mat_vasc = vasc{i}; % Snapshot from red channel

    % Iterate through the rows in the feature matrix
    for j=1:size(mat,1)
        
        % Define classes for morphological classification
        classNames = {'Branching','Diffuse translocation', 'Junk','Locomotion', 'Perivascular translocation', 'Round'};
        
        % Create an RGB image from the snapshots
        im = zeros([blocksize blocksize 3]);
        im(:,:,2) = reshape(mat(j,:), [blocksize blocksize]);
        im(:,:,1) = reshape(mat_vasc(j,:), [blocksize blocksize]);
        
        % Perform channel-wise mean subtraction centering
        im_predict = im/255 - mean_image_6class;
        im_predict_tme = im/255 - mean_image_tme;

        % Predict morphology
        morph_class = morph_net.predict(im_predict);
      
        % Class with highest probability is chosen
        [~, predictedIndex_61] = max(morph_class);
        predictedIndex = predictedIndex_61;
        predictedLabel = onehotdecode(morph_class,classNames,2);

        % In debug mode snapshots are written to a folder
        if(debug)
            if(max(morph_class) < 0.6)
                folder = [char(predictedLabel) '/<60%/'];
                filename = ['score_' num2str(max(morph_class)) '_stack_' num2str(m) '_im_' num2str(i) '_' num2str(randi([1,1000])) '_' '.tif'];
                imwrite(uint8(im), [output_path folder filename],'tiff')
            else
                folder = [char(predictedLabel) '/>60%/'];
                filename = ['score_' num2str(max(morph_class)) '_stack_' num2str(m) '_im_' num2str(i) '_' num2str(randi([1,1000])) '_' '.tif'];
                imwrite(uint8(im), [output_path folder filename],'tiff')
            end
        end
       
        % Save class probabilities and label
        regs = array2table(morph_class, 'VariableNames',classNames);
        regs.predicted_label = predictedIndex;

        % Define TME classes
        classNames2 = {'Microglia colocalized', 'Vessel associated'};

        % Predict TME interaction
        tme_interact = trainedNetwork_tme.predict(im_predict_tme);
        [score,predictedIndex_tme] = max(tme_interact);
        predictedLabel_tme = onehotdecode(tme_interact,classNames2,2);

        % Threshold the score so that if lower than some thresholds,
        % the image is considered to represent no interaction
        if(predictedLabel_tme == "Microglia colocalized" && score < 0.975)
            folder = 'Non-associated/';
            filename = ['score_' num2str(max(tme_interact)) '_stack_' num2str(m) '_im_' num2str(i) '_' num2str(randi([1,1000])) '_' '.tif'];
            predictedIndex_tme = 3;
        elseif(predictedLabel_tme == "Vessel associated" && score < 0.8)
            folder = 'Non-associated/';
            filename = ['score_' num2str(max(tme_interact)) '_stack_' num2str(m) '_im_' num2str(i) '_' num2str(randi([1,1000])) '_' '.tif'];
            predictedIndex_tme = 3;
        else
            folder = [char(predictedLabel_tme) '/'];
            filename = ['score_' num2str(max(tme_interact)) '_stack_' num2str(m) '_im_' num2str(i) '_' num2str(randi([1,1000])) '_' '.tif'];
        end

        % Save TME label
        regs.TME_status = predictedIndex_tme;
        
        % If in debug mode, display image + label and write snapshot to folder
        if(debug)
            imagesc(im/255)
            title(char(predictedLabel));
            imwrite(uint8(im), [output_path folder filename],'tiff')
        end
        
        % Save predicted features to stucture
        properties{i}(j,:) = regs;
    end
    
end

end