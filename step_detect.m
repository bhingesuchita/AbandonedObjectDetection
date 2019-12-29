function [comp,frame] = step_detect(time_signatures,varargin)
% Function to detect abandoned objects. The time signature of an abandoned
% object has a step-type response. 
% Input:
%   time_signatures : matrix of dimension nFrames X order (M_bar) 
%   c1 : correlation threshold in stage 1 to identify time points of potential step
%   change (default : 0.8)
%   c2 : correlation threshold in stage 2 to identify step-type response
%   (default : 0.8)
% Output:
%   comp : cell-structure containing component index containing step-type response
%   frame : cell-structure containing frame index at which step increase occurred

c1 = 0.8;c2 = 0.8; % set default parameters

% Read user-defined parameters
if mod(length(varargin),2) == 1
    error("Input arguments must be in name-value pairs");
else
    i = 1;
    while i <= length(varargin)
        if varargin{i} == "c1"
            c1 = varargin{i+1};
        elseif varargin{i} == "c2"
            c2 = varargin{i+1};
        end
        i = i + 2;
    end
end

% Get dimensions of input matrix
[nFrames,order,K] = size(time_signatures);

% Define ideal step of length 200 for correlation
ffilter = [-1*ones(1,100) ones(1,100)];

for k = 1 : K % loop through all datasets
    z = zscore(time_signatures(:,:,k)); % normalize
    
    % Initialize variables for later use
    comp_k = []; 
    frame_k = []; 
    for m = 1 : order % loop through all components
        
        % Stage 1 : Identify points of potential step changes/ regions of interest
        tvals = 0; 
        for j = 1 : nFrames-200
            y(j) = corr2(ffilter.',z(j:j+199,m));
            if((abs(corr2(ffilter.',z(j:j+199,m))) > c1)) % Select regions-of-interest based on threshold c1
                [~,~,~,stats] = ttest2(z(1:j+100,m),z(j+101:nFrames,m));
                if(abs(stats.tstat) > tvals)
                    index = j+100; % store index corresponding to highest t-value
                    tvals = abs(stats.tstat); % store maximum t-value
                end
            end     
        end
        
        % Stage 2 : In order to identify if the time signature is a step-type
        % response, we perform k-means clustering
        if(tvals ~= 0)
            [label,~] = kmeans(z(:,m),2);
            f = [ones(1,index) 2*ones(1,nFrames-index)];
            if(abs(corr2(label.',f)) > c2) %c2 = .8
                comp_k = [comp_k m]; % store the index of the step-type response
                frame_k = [frame_k index]; % store the frame index at which the step change occurred
            end
        end
    end
    comp{k} = comp_k;
    frame{k} = frame_k;
end
