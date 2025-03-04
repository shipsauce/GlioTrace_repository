function combine_video_subroutine(video_files,output_filename,frame_rate,output_width,output_height,legends,slowdown,videoquality,gammavalues)
% This function takes a set of videos and combines them into one
% video, which is written to an output directory.
%
% Input parameters:
%  video_files - paths to individual videos
%  output_filename - name of combined video
%  frame_rate - frames per second
%  output_width - width of output video
%  output_height - height of output video
%  legends - titles and subtitles of idnvidual videos
%  slowdown - vector of values for each video, whether to slowdown the
%                  frame-rate
%  videoquality - resolution of output video
%  gammavalues - values for gamma correction
%
% @authors: Sven Nelander
% @date: 28022025
%

% Define output video parameters

%output_width = 1920; % Standard 16:9 width
%output_height = 1080; % Standard 16:9 height

text_position = [50, output_height*0.83]; % Position (X, Y) in pixels
text_position2 = [50, output_height*0.9]; % Position (X, Y) in pixels

text_color = [255, 255, 255]; % White color
text_font_size = 30; % Font size

num_videos = length(video_files);

% Create video writer object
output_video = VideoWriter(output_filename, 'MPEG-4');
output_video.FrameRate = frame_rate;
   output_video.Quality = videoquality;
open(output_video);

% Process each video
for i = 1:num_videos
    'video' 
    i
    % Read video
    video_reader = VideoReader(video_files{i});
    input_width = video_reader.Width;
    input_height = video_reader.Height;

    % Calculate aspect ratios
    input_aspect_ratio = input_width / input_height;
    target_aspect_ratio = output_width / output_height;

    % Determine scaling factor
    if input_aspect_ratio > target_aspect_ratio
        % Wider than 16:9 - Fit width, adjust height
        new_width = output_width;
        new_height = round(output_width / input_aspect_ratio);
    else
        % Taller than 16:9 - Fit height, adjust width
        new_height = output_height;
        new_width = round(output_height * input_aspect_ratio);
    end

    % Positioning on black background
    x_offset = round((output_width - new_width) / 2);
    y_offset = round((output_height - new_height) / 2);

    % Create output frames
    while hasFrame(video_reader)
        % Read frame
        frame = readFrame(video_reader);
        
        green=frame(:,:,2);
        green=(double(green)/255).^gammavalues(i);
        green=uint8(green*255);
        frame(:,:,2)=green;

        % Resize while maintaining aspect ratio
        resized_frame = imresize(frame, [new_height, new_width]);

        % Create black background
        black_frame = zeros(output_height, output_width, 3, 'uint8');

        % Place resized video at the center
        black_frame(y_offset+1:y_offset+new_height, x_offset+1:x_offset+new_width, :) = resized_frame;

                black_frame = insertText(black_frame, text_position, legends{i}{1}, ...
                                 'FontSize', 1.5*text_font_size, ...
                                 'BoxColor', 'black', 'TextColor', text_color, ...
                                 'BoxOpacity', 0);
                    black_frame = insertText(black_frame, text_position2, legends{i}{2}, ...
                                 'FontSize', text_font_size, ...
                                 'BoxColor', 'black', 'TextColor', text_color, ...
                                 'BoxOpacity', 0);
        % Write frame to output video
        for reps=1:slowdown(i)
            writeVideo(output_video, black_frame);
        end

    end

    for reps=1:19
  
        writeVideo(output_video, uint8(double(black_frame)*(20-reps)/20));
    end

end

% Close video writer
close(output_video);

disp('Combined video saved successfully!');




%%


