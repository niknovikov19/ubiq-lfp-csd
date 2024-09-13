%% plots alpha aligned csd and mua

clear
close all

% make this 1 if you want to see acoustic responses 
gimmeplots=1;

% select trigger, 1 for aud, 3 for visual, 4 for sacc on, 5 for sacc off
trigch=1;

%% path info 

paths = {'E:\spectrolaminar\spont\core\1-kk059060048@os.mat'
};

%'E:\spectrolaminar\spont\core\1-kk041042049@os.mat' 'E:\spectrolaminar\spont\core\1-kk059060048@os.mat' 
% Load data
load(paths{1});  % Assume craw is loaded

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

% Extract continuous ephys data and triggers
[~, cnte, cntm, cntc, ~, cntb] = module_cnt05(craw, newadrate, filtere, filteru, filtertype);

%% Bandpass filter for CSD (8-14 Hz)

[b, a] = butter(4, [8 14] / (newadrate / 2), 'bandpass'); %  Butterworth filter

% Apply bandpass filter
bpCntc = filtfilt(b, a, cntc')'; 

% Find troughs in bandpass filtered data
centralChannelData = bpCntc(16, :);
centralChannelStd = std(centralChannelData)*2;
[~, locs] = findpeaks(-centralChannelData, 'MinPeakHeight', centralChannelStd);

if length(locs)<100
    disp('less than 100 troughs found!')
end
%% Epoch data according to epoch_tframe
% Initialize variables
csdEpochs = [];
muaEpochs = [];
nChannels = size(cntc, 1);
nTrials = length(locs);

for i = 1:nTrials
    startIdx = locs(i) + epoch_tframe(1) * newadrate / 1000;
    endIdx = locs(i) + epoch_tframe(2) * newadrate / 1000;

    if startIdx > 0 && endIdx <= size(cntc, 2)
        csdEpochs(:, i, :) = cntc(:, startIdx:endIdx);
        muaEpochs(:, i, :) = cntm(:, startIdx:endIdx);
    end
end

%% Remove outliers
[csdEpochsClean, muaEpochsClean,~] = rejectartifacts_CnMonly(csdEpochs,muaEpochs);
% csdEpochsClean = csdEpochs;
% muaEpochsClean = muaEpochs;
%% Plotting if required
if gimmeplots
    figure;
    % Average CSD and MUA across trials
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
    caxis([5.25 7.25])
end
