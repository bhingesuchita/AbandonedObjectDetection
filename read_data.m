function [dataOut,m,n] = read_data(folder,skip_frames)

%%%%%%%%%%%%%% Code description
% Input:
%   folder : folder containing input data in the form of sequence of images
%   ('.jpg','.png','.tif') or video files
%           
%       Video from multiple cameras:
%           If data is in the form of sequence of images, store them in
%           different subfolders
%           Output data will be three dimensional of the form T X P X K,
%           where T is the number of frames, P is the number of pixels and
%           K is the number of cameras/video sequences
%   
%   skip_frame : frame rate. Reads every 5th (default) frame. Set 1 to read
%   all frames
% 
% Output:
%   dataOut : video/image data converted to gray-scaled images and stored in a matrix.
%       Video from single camera:
%           Output data will be two dimensional of the form T/skip_frames X P,
%           where T is the number of frames and P is the number of pixels 
%       Video from multiple cameras:
%           Output data will be three dimensional of the form T/skip_frames X P X K,
%           where T is the number of frames, P is the number of pixels and
%           K is the number of cameras/video sequences. 
%   m,n : size of image in the video
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imds = imageDatastore(folder,'IncludeSubfolders',1,'LabelSource','foldernames'); % form image datastore
[m,n,~] = size(imread(imds.Files{1}));
imdsgray = transform(imds,@(x) transform_RGB2gray(x,skip_frames)); % transform all RGB images to gray-scale

countLabel = countEachLabel(imdsgray.UnderlyingDatastore); % count number of frames in each subfolder, T should be same for all K
K = size(countLabel,1); % number of subfolders, i.e., number of cameras

if sum(countLabel.Count./countLabel.Count(1)) ~= K % check if number of frames in each subfolder are equal
    error("Number of frames from all cameras must be equal");
else
    T = countLabel.Count(1);
end

imdsgray.UnderlyingDatastore.ReadSize = T; % specify the number of files to read equal to number of frames
clear imds
if K == 1 % Video files from a single camera
    data = read(imdsgray);
elseif K > 1 % Video files from multiple cameras
    k = 1; % counts the number of subfolders/cameras
    while hasdata(imdsgray)
        data = read(imdsgray); % reads top T images from datastore and removes them
        
        % reshape each image to form a vector
        dataOut(:,:,k) = cell2mat(transform_2DImage2vector(data,m,n));
        k = k + 1; % increment subfolder counter
    end
end
end

function dataOut = transform_RGB2gray(dataIn,skip_frames)
% Input : 
%   dataIn : RGB image data stored in cell format of size T X 1. Each cell
%   contains image data of size m x n x 3, where m and n are dimensions of
%   the image.
%   skip_frames : frame rate. Selects every 5th (default) frame
%
% Output:
%   dataOut : gray-scaled images stored in cell format of size
%   T/skip_frames X 1. Each cell contains image data of size m X n.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = 1;
    for i = 1 : skip_frames : length(dataIn)
        dataOut{j,1} = rgb2gray(dataIn{i,1});
        j = j + 1;
    end
end

function [dataOut] = transform_2DImage2vector(dataIn,m,n)
% Input:
%   dataIn : gray-scaled images stored in cell format of size
%   T/skip_frames X 1. Each cell contains image data of size m X n.
%
% Output:
%   dataOut : vectorised images stored in cell format. Each cell contains
%   image data of size 1 x m*n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(dataIn)
    dataOut{i,1} = reshape(dataIn{i,1},1,m*n);
end
end