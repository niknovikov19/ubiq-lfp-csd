%% FILTERED LFP
clear
clc
close all
load data.mat;
%%
fs = 1000; % Sampling rate in Hz

% Butterworth low-pass filter from 0.1 to 300 Hz
[b, a] = butter(2, [0.1 300] / (fs / 2));  % Use 'low' by default, as this is a bandpass

% Select the LFP data for example1
%lfpData = data.example1_vlPFC_lfp; % channels x trials x timePoints
% path = 'E:\spectrolaminar\AttnData\core\cont\1-ab025026023@os.mat';
path = 'E:\spectrolaminar\AttnData\core\cont\1-kk041042029@os.mat';
[~, ~, ~, ~, ~,~,~,lfpData] = EphysExtractFxn(path,1);


numChannels = size(lfpData, 1);
numTrials = size(lfpData, 3);
numTimePoints = size(lfpData, 2);
% Time vector assuming a sampling rate of 1,000 Hz
timeVector = linspace(0, (numTimePoints - 1) / fs, numTimePoints);

% Initialize matrix for filtered data (channel x trial x time)
% filteredLFP = nan(numChannels, numTrials, numTimePoints);
filteredLFP = lfpData;

% % Filter data
% for channel = 1:numChannels
%     for trial = 1:numTrials
%         tempTrialData = squeeze(lfpData(channel, trial, :));
%         if ~all(isnan(tempTrialData))
%             % Apply Butterworth low-pass filter
%             filteredLFP(channel, trial, :) = filtfilt(b, a, tempTrialData);
%             
%         end
%     end
% end

%% csd
CSD=[];
CSD(:,:,:) = -diff(filteredLFP,2,1); %CSD


%% Next nearest neighbor CSD calculation
CSDnn = nan(numChannels-4, numTrials, numTimePoints); % Allocate space for CSD
for trial = 1:numTrials
    for t = 1:numTimePoints
        for ch = 3:numChannels-2
            CSDnn(ch-2, trial, t) = (filteredLFP(ch+2, trial, t) - 2 * filteredLFP(ch, trial, t) + filteredLFP(ch-2, trial, t)) / (4 * (1/fs)^2);
        end
    end
end

%% wavelet

Fs = 1000;
frequencyRange = [1 40];
targetBand = [10 30];
% Initialize storage for results
% Initialize storage for amplitude and phase for each channel, trial, and time
powerLFP = zeros(numChannels, numTrials, numTimePoints);
phaseLFP = zeros(numChannels, numTrials, numTimePoints);

powerCSD = zeros(numChannels-2, numTrials, numTimePoints);
phaseCSD = zeros(numChannels-2, numTrials, numTimePoints);

% Perform wavelet analysis on each trial of each channel
for channel = 1:numChannels
    for trial = 1:numTrials
        % Perform analysis on LFP
        [amplitudeLFP, phaseLFPData] = wvlt_bndlm_fxn(squeeze(filteredLFP(channel, :, trial)), Fs, frequencyRange,targetBand);
        powerLFP(channel, trial, :) = amplitudeLFP;
        phaseLFP(channel, trial, :) = phaseLFPData;

        % Perform analysis on CSD if the channel is within bounds
        if channel <= numChannels - 2
            [amplitudeCSD, phaseCSDData] = wvlt_bndlm_fxn(squeeze(CSD(channel, :, trial)), Fs, frequencyRange,targetBand);
            powerCSD(channel, trial, :) = amplitudeCSD;
            phaseCSD(channel, trial, :) = phaseCSDData;
        end
    end
end


% Assuming powerLFP, phaseLFP, powerCSD, and phaseCSD are already computed and stored

%% fig


% Assuming filteredLFP and CSD have been computed and stored

%% Create Heatmaps for LFP Data
figure; % Create a new figure window for LFP
subplot(3,2,1)
imagesc(timeVector, 1:numChannels, squeeze(mean(filteredLFP, 3)));
title('Trial Avg. LFP');
xlabel('Time (s)');
ylabel('Channel');
colorbar; % Show a color bar to indicate the scale
colormap jet; % Set colormap
%caxis([min(filteredLFP(:))*0.75, max(filteredLFP(:))*0.1]); % Adjust color axis for better visualization

% Create Heatmaps for CSD Data % Create a new figure window for CSD
subplot(3,2,2)
imagesc(timeVector, 1:numChannels-2, squeeze(mean(CSD, 3)));
title('Trial Avg. CSD');
xlabel('Time (s)');
ylabel('Channel');
colorbar;
colormap jet;
%caxis([min(CSD(:)), max(CSD(:))]); % Adjust color axis for better visualization
%
 % Create a new figure window

% Heatmap for LFP Power
subplot(3, 2, 3); % Positioning this plot in a 2x2 grid
imagesc(timeVector, 1:numChannels, squeeze(mean(powerLFP, 3)));
title('LFP Alpha-Beta Abs. Amp');
xlabel('Time (s)');
ylabel('Channel');
colorbar; % Show a color bar to indicate the scale
colormap jet; % Set colormap

% Heatmap for LFP Phase
subplot(3, 2, 5);
imagesc(timeVector, 1:numChannels, squeeze(mean(phaseLFP, 3)));
title('LFP Alpha-Beta Amp.');
xlabel('Time (s)');
ylabel('Channel');
colorbar;
colormap jet;

% Create Heatmaps for CSD Data
 % Create a new figure window

% Heatmap for CSD Power
subplot(3, 2, 4);
imagesc(timeVector, 1:numChannels-2, squeeze(mean(powerCSD, 3)));
title('CSD Alpha-Beta Abs. Amp');
xlabel('Time (s)');
ylabel('Channel');
colorbar;
colormap jet;

% Heatmap for CSD Phase
subplot(3, 2, 6);
imagesc(timeVector, 1:numChannels-2, squeeze(mean(phaseCSD, 3)));
title('CSD Alpha-Beta Amp.');
xlabel('Time (s)');
ylabel('Channel');
colorbar;
colormap jet;

%% Compute and plot averaged power over time as a function of channel

% Average power over time and trials for LFP
% avgPowerLFP = squeeze(mean(mean(powerLFP, 3), 2));  % Average across time (3rd dim) and trials (2nd dim)
% 
% % Average power over time and trials for CSD
% avgPowerCSD = squeeze(mean(mean(powerCSD, 3), 2));  % Same, adjust for reduced channel count in CSD
% 
% %
% normLFP = avgPowerLFP/max(avgPowerLFP);
% normCSD = avgPowerCSD/max(avgPowerCSD);
% normLFPg = avgGammaPowerLFP/max(avgGammaPowerLFP);
% normCSDg = avgGammaPowerCSD/max(avgGammaPowerCSD);
% 
% % Create a new figure for these plots
% figure;
% subplot(2, 1, 1);  % LFP Power Plot
% plot(normLFP,1:numChannels,  'LineWidth', 2);
% hold on
% plot(normLFPg,1:numChannels,  'LineWidth', 2);
% set(gca, 'YDir','reverse');  % Reverse the direction of the y-axis
% legend('AB','Gamma')
% title('Average Power LFP ');
% ylabel('Channel');
% xlabel('Average Power');
% grid on;
% 
% 
% % Enhance layout
% set(gca, 'FontSize', 12);  % Set font size for readability
% 
