%% Initialization
clear;
close all;

%% user inputs
%baselinePeriod = [-200 -100]; % Baseline in ms, comment out to inactivate
epoch_tframe = [-250 250];  
newadrate = 1000;
filtere = [0.5 300]; % LFP
filteru = [300 5000]; % MUA

%%
% Directory and parameters
dataDir = 'E:\spectrolaminar\spont\core\'; % Directory containing .mat files
figuresDir = fullfile(dataDir, 'muaChans'); % Directory to save results
if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

% Get all files with "@os" in their name
filePattern = '*@os*.mat';
myfiles = dir(fullfile(dataDir, filePattern));
% Initialize the structure to collect all results
allResults = struct('fileName', {}, 'channelMeasureIdx', {});
%% Loop through each file
for fileIdx = 1:length(myfiles)
    % Load the .mat file
    basefilename = myfiles(fileIdx).name;
    fullfilename = fullfile(dataDir, basefilename);
    load(fullfilename); % Assume it loads `craw`
    
    % Define where to save the data
    saveFilename = fullfile(figuresDir, [basefilename(1:end-4), '_results.mat']);
    
    %% INPUTS (filtering, epoch timing etc.)

    filtertype = 1;
    fsize = 6;
    xlabelres = 5;

    %% Set timeframe
    time = epoch_tframe(1):1000/newadrate:epoch_tframe(2);
    zeroIndex = find(time == 0); % Find index for 0 ms (midpoint of epoch)
    % Determine baseline period indices
    if exist('baselinePeriod','var')
        baselineIndices = find(time >= baselinePeriod(1) & time <= baselinePeriod(2));
    end

    % Extract continuous ephys data and triggers
    [~, cnte, cntm, cntc, ~, cntb] = module_cnt05(craw, newadrate, filtere, filteru, filtertype);

    %% Bandpass filter for CSD (8-14 Hz)
    [b, a] = butter(4, [8 14] / (newadrate / 2), 'bandpass'); % Butterworth filter
    bpCntc = filtfilt(b, a, cntc')'; 

    %% Initialize storage for MUA amplitudes at 0 ms
    numChannels = size(cntc, 1);
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
                if exist('baselinePeriod','var')
                    baselineMean = mean(muaEpochsClean(ch, :, baselineIndices), 'all');
                    muaAmplitudes(channel, ch) = mean(muaEpochsClean(ch, :, zeroIndex), 'all')-baselineMean;
                else
                    muaAmplitudes(channel, ch) = mean(muaEpochsClean(ch, :, zeroIndex), 'all');
                end
                    
            end
        end
    end

    %% Find the channel combination with the maximum MUA amplitude
    [~, linearIndex] = max(muaAmplitudes(:));
    [channelAlignIdx, channelMeasureIdx] = ind2sub(size(muaAmplitudes), linearIndex);

    %% Store results
    saveData.fileName = basefilename;  % Store the file name
    saveData.channelMeasureIdx = channelMeasureIdx;

    %% Save the data
    save(saveFilename, '-struct', 'saveData');
    
    %% Collect results for final summary
    allResults(fileIdx).fileName = basefilename;
    allResults(fileIdx).channelMeasureIdx = channelMeasureIdx;

    % Display progress
    fprintf('Processed and saved: %s\n', basefilename);
end
%% Save all results to a single file
finalSaveFilename = fullfile(figuresDir, 'all_file_results.mat');
save(finalSaveFilename, 'allResults');

fprintf('All results saved to: %s\n', finalSaveFilename);