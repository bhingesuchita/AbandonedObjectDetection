function detect_abandonedObject(folder,varargin)

%%%%%%%%%%%%%% Code description
% Input:
%   folder : folder containing input data in the form of sequence of images
%   ('.jpg','.png','.tif') or video files
%
% varargin - an optional input which is a structure containing the
% following elements. Note that default values are listed next to the 
% structure elements.
%
%    'skip_frames',5, ... % frame rate. Reads every 5th (default) frame.
%    Set 1 to read all frames
%    'nRuns',10, ... % Number of runs for ICA/IVA

% Please cite following references:
%   [1] Y.-O. Li, T. Adali, and V. D. Calhoun, “Estimating the number of 
%       independent components for functional magnetic resonance imaging 
%       data,” Human brain mapping, vol. 28, no. 11, pp. 1251–1266, 2007.
%   [2] Xi-Lin Li, Tulay Adali, "Blind spatiotemporal separation of second
%       and/or higher-order correlated sources by entropy rate 
%       minimization," IEEE International Conference on Acoustics, Speech
%       and Signal Processing 2010.
%   [3] M. Anderson, G.-S. Fu, R. Phlypo, and T. Adali, “Independent vector
%       analysis, the Kotz distribution, and performance bounds,” in IEEE 
%       International Conference on Acoustics, Speech and Signal Processing
%       (ICASSP), 2013, pp. 3243–3247.
%   [4] Q. Long, C. Jia, Z. Boukouvalas, B. Gabrielson, D. Emge and T.
%       Adali, "Consistent Run Selection for Independent Component Analysis:
%       Application to FMRI Analysis," in International Conference on 
%       Acoustics, Speech and Signal Processing (ICASSP), Calgary, AB, 
%       2018, pp. 2581-2585.
%   [5] S. Bhinge, Y. Levin-schwartz, and T. Adali, "Data-driven fusion of
%       multi-camera video sequences: Application to abandoned object
%       detection," in IEEE International Conference on Acoustics, Speech
%       and Signal Processing (ICASSP), 2017, pp. 1697-1701.
%
% Code written by : Suchita Bhinge (suchita1@umbc.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gather default options structure
options=struct( ...
   'skip_frames',5, ... % frame rate. Reads every 5th (default) frame. Set 1 to read all frames
   'nRuns',2 ... % Number of runs for ICA/IVA
   );
% load in user supplied options
options=getopt(options,varargin{:});

% Read data
[data,m,n] = read_data(folder,options.skip_frames); % data is of dimension
data = double(data);
% nFrames X nPixels X nCameras, m and n are the dimensions of the images

% Get dimension of data
[~, ~, nCameras] = size(data); 

% Model order selection : 
% The problem is overdetermined in nature, i.e, the
% dimension of the signal subspace is less than the number of
% frames/observation. Hence, we perform dimension reduction to obtain the
% signal subspace and the mode order, dimension of signal subspace, is 
% estimated using E-DS method [1]

for k = 1 : nCameras
    order(k) = OrderSelection_E_DS(data(:,:,k)');
end

% Get final order, M_bar, as median across orders from all cameras
M_bar = median(order);

% Perform PCA to estimate the signal subspace for each camera
if nCameras == 1
    [coeffX,temp] = pca(data(:,:,k)','NumComponents',M_bar); 
    % The columns of temp are the principal components and coeffX is the
    % data reduction matrix
    
    X = temp'; % The rows of X are the principal components
else
    for k = 1 : nCameras
        [coeffX(:,:,k),temp] = pca(data(:,:,k),'NumComponents',M_bar); 
        % The columns of temp are the principal components and coeffX is the
        % data reduction matrix
    
        X(:,:,k) = temp'; % The rows of X are the principal components 
    end
end

% Perform ICA using ERBM algorithm [2] if nCameras = 1 and transposed-IVA
% using IVA-GGD algorithm [3] if nCameras > 1.
% The code for ERBM and IVA-GGD algorithms can be obtained from 
% "http://mlsp.umbc.edu/resources.html" We typically perform multiple 
% ICA/IVA runs with different initial points and select the consistent run
% using [4]. 

if nCameras == 1
    for run = 1 : options.nRuns
        fprintf('Run %d\t',run);
        W{run} = ERBM(X,11); % Perform ERBM with filter length set to 11 (default)
        fprintf('complete\n');
    end
else
    for run = 1 : options.nRuns
        fprintf('Run %d\t',run);
        W{run} = iva_mpe_decp_v2(X,'verbose',true,'alpha0',2.0); % Perform IVA-GGD (also known as IVA-MPE-decp)
        fprintf('complete\n');
    end
end

% Consistent run selection using [4]
if nCameras == 1
    selRun = RunSelection_crossISIidx(W);
    W_best = W{selRun};
else
    selRun = run_selection_crossJointISI(W);
    W_best = W{selRun};
end

% Compute the mixing matrices and sources
for k = 1 : nCameras
    A(:,:,k) = coeffX(:,:,k)*inv(W_best(:,:,k)); 
    S(:,:,k) = W_best(:,:,k) * X(:,:,k); 
end

% If nCameras = 1 the columns of A are the time signatures and rows of S
% are the sources representing objects in the video. If nCameras > 1, we
% use the transposed model. Hence the columns of A are the sources and rows
% of S are the time signatures. We use the time signatures to detect the
% abandoned object.
if nCameras == 1
    [comp,frame_index] = step_detect(A);
else
    [comp,frame_index] = step_detect(S);
end

% Plot detected objects and estimated frame index at which object was
% dropped
if nCameras == 1
    for i = 1 : length(comp)
        S(comp(i),:) = (S(comp(i),:)-min(S(comp(i),:)))/(max(S(comp(i),:)) - min(S(comp(i),:))); % Normalize between 0 and 1 for plotting
        figure();
        imshow(reshape(S(comp,:),m,n)); % Reshape and plot
        title(sprintf('Abandoned object detected from Camera %d at frame %d',i,frame_index(i)));
    end
else
    for k = 1 : nCameras
        for i = 1 : length(comp{k})
            A(:,comp{k}(i),k) = (A(:,comp{k}(i),k)-min(A(:,comp{k}(i),k)))/(max(A(:,comp{k}(i),k)) - min(A(:,comp{k}(i),k))); % Normalize between 0 and 1 for plotting
            figure();
            imshow(reshape(A(:,comp{k}(i),k),m,n));
            title(sprintf('Abandoned object detected from Camera %d at frame %d',k,frame_index{k}(i)));
        end
    end
end