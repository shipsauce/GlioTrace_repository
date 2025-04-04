
clc; clear; close all;

% -----------------
% User Settings
% -----------------
user_specific_path='define_path';
folderpath = [user_specific_path 'Supplementary_videos_MP4/']; % Change to your folder path
output_filename = [user_specific_path 'Supplementary_videos_MP4/comb.mp4']; % Change output path
titles = {"Supplementary video 1", "TL-DyLight-594 bioavailability"; "Supp 2", "Description"; "Supp 3", "Description"}; % Titles corresponding to each video
fps = 7; % Frame rate of the output video
titleDisplayTime = 2.5; % Time in seconds the title is shown before the video starts

% -----------------
% Read Video Files
% -----------------
videoFiles = dir(fullfile(folderpath, '*.mp4'));
nVideos = length(videoFiles);

if nVideos == 0
    error('No MP4 files found in the folder.');
end


%%

% Name of output video files
final_video_files={'Supplementary_Video_1A-1G.mp4','Supplementary_Video_2A-2F.mp4','Supplementary_Video_3A-3G.mp4','Supplementary_Video_4A-4E.mp4'};

% Load metadata table
t=readtable('videolist.xlsx');
video_quality_tweak=[75 100 74 70];

% Create the output videos
for iteration=1:4

    % Merge output folder and output video title
    output_filename = [user_specific_path final_video_files{iteration}];

    % Find the rows of the original videos that should be included in the merged
    % video
    f=find(t.video_index==iteration);

    % Create paths to the origial videos
    input_videos={};
    for j=1:length(f)
        input_videos{j,1}=[folderpath t.original_file{f(j)}];
    end
        
    % Find metadata for original videos
    legends={};
    for j=1:length(f)
        legends{j}={t.label_in_paper{f(j)} , t.description{f(j)}}
    end

    % Define video output parameters
    slowdown=t.slowdown(f);
    gammavalues=t.gamma(f);
    
    % Create merged video
    combine_video_subroutine(input_videos,output_filename,24,1920,1080,legends,slowdown,video_quality_tweak(iteration),gammavalues)

end