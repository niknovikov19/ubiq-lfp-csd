%% plots alpha aligned csd and mua

clear
close all

% make this 1 if you want to see acoustic responses 
gimmeplots=0;

% select trigger, 1 for aud, 3 for visual, 4 for sacc on, 5 for sacc off
trigch=1;

%% path info 

paths = {'E:\spectrolaminar\spont\core\1-kk059060048@os.mat'};
% Load data
load(paths{1});  % Assume craw is loaded
% Define directory to save the data
figuresDir = 'E:\spectrolaminar\spont\core\muaChans\'; 


%% INPUTS (filtering, epoch timing etc.)%

epoch_tframe = [-250 250];  
newadrate = 1000;
filtere = [0.5 300]; % LFP
filteru = [300 5000]; % MUA
filtertype = 1;
fsize = 6;
xlabelres = 5;

%% Set timeframe
time = epoch_tframe(1):1000/newadrate:epoch_tframe(2);
zeroIndex = find(time == 0); % Find index for 0 ms (midpoint of epoch)

% Extract continuous ephys data and triggers
[~, cnte, cntm, cntc, ~, cntb] = module_cnt05(craw, newadrate, filtere, filteru, filtertype);

%% Bandpass filter for CSD (8-14 Hz)
[b, a] = butter(4, [8 14] / (newadrate / 2), 'bandpass'); % Butterworth filter

% Apply bandpass filter
bpCntc = filtfilt(b, a, cntc')'; 

%% Initialize storage for MUA amplitudes at 0 ms
numChannels = size(cntc,1);
muaAmplitudes = zeros(numChannels, numChannels); % Preallocate a matrix to store MUA amplitudes for each channel alignment

%% Loop through each channel for trough detection and MUA alignment
for channel = 1:numChannels
    % Extract data for the current channel
    currentChannelData = bpCntc(channel, :);
    
    % Compute the standard deviation for the current channel
    currentChannelStd = std(currentChannelData)*2;
    
    % Find troughs in the current channel's bandpass-filtered data
    [~, locs] = findpeaks(-currentChannelData, 'MinPeakHeight', currentChannelStd);
    
    %% Epoch data according to epoch_tframe
    csdEpochs = [];
    muaEpochs = [];
    nChannels = size(cntc, 1);
    nTrials = length(locs);
    
    for i = 1:nTrials
        startIdx = locs(i) + round(epoch_tframe(1) * newadrate / 1000);
        endIdx = locs(i) + round(epoch_tframe(2) * newadrate / 1000);

        if startIdx > 0 && endIdx <= size(cntc, 2)
            csdEpochs(:, i, :) = cntc(:, startIdx:endIdx);
            muaEpochs(:, i, :) = cntm(:, startIdx:endIdx);
        end
    end
    
    %% Remove outliers
    [csdEpochsClean, muaEpochsClean] = rejectartifacts_CnMonly(csdEpochs, muaEpochs);

    %% Extract MUA amplitude at 0 ms for all channels
    if ~isempty(muaEpochsClean)
        % Mean MUA amplitude across trials at 0 ms
        for ch = 1:numChannels
            muaAmplitudes(channel, ch) = mean(muaEpochsClean(ch, :, zeroIndex), 'all');
        end
    end
end


% Find the maximum MUA amplitude and its indices
[maxMUAamp, linearIndex] = max(muaAmplitudes(:));

% Convert linear index to row and column indices
[channelAlignIdx, channelMeasureIdx] = ind2sub(size(muaAmplitudes), linearIndex);

% Store path and channelMeasureIdx
saveData.path = paths{1};  % Use paths{1} to save the first path
saveData.channelMeasureIdx = channelMeasureIdx;

% Save the data
save(fullfile(figuresDir, 'path_and_channelMeasureIdx.mat'), '-struct', 'saveData');

%% Plotting if required
if gimmeplots
    figure;
    % Average CSD and MUA across trials for the last channel used
    meanCSD = mean(csdEpochsClean, 2);
    meanMUA = mean(muaEpochsClean, 2);
    
    % Plot CSD
    subplot(2,1,1);
    imagesc(time, 1:nChannels, squeeze(meanCSD));
    title('Average CSD');
    xlabel('Time (ms)');
    ylabel('Channel');
    colorbar;
    
    % Plot MUA
    subplot(2,1,2);
    imagesc(time, 1:nChannels, squeeze(meanMUA));
    title('Average MUA');
    xlabel('Time (ms)');
    ylabel('Channel');
    colorbar;
end
