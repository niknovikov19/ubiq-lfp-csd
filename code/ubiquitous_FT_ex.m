% find spectrolaminar trends using field trip
% chase m 2024

clear
close all

% data3D with dimensions channels, trials, epoched time
load data.mat
lfp1 = data.example1_vlPFC_lfp;
% path = 'E:\spectrolaminar\AttnData\core\cont\1-ab025026023@os.mat';
% [~, ~, ~, ~, ~,~,~,lfp1] = EphysExtractFxn(path,1);
data3D = lfp1;
data = struct();
data.label = {};      % Update this with your channel names
data.fsample = 1000;  %  sampling frequency

% Number of channels, trials, and time points
[numChannels, numTimePoints, numTrials] = size(data3D);

% Set channel labels (assuming generic labels here)
data.label = arrayfun(@(x) sprintf('Ch%d', x), 1:numChannels, 'UniformOutput', false);

% Reshape 3D double array data for field trip
data.trial = cell(1, numTrials);
data.time = cell(1, numTrials);

% Populate the trial and time fields
for trial = 1:numTrials
    % Each cell in data.trial should be [numChannels x numTimePoints]
    % Note: You might need to transpose this depending on your data orientation
    data.trial{trial} = squeeze(data3D(:, :, trial));
    
    % Assuming continuous and uniform time sampling across trials
    data.time{trial} = (0:numTimePoints-1) / data.fsample;
end

% Want to only analyze part of the epoch?
numTimePointsToAnalyze = numTimePoints; % Number of time points to analyze
for trial = 1:length(data.trial)
    if length(data.time{trial}) >= numTimePointsToAnalyze
        data.trial{trial} = data.trial{trial}(:, 1:numTimePointsToAnalyze);
        data.time{trial} = data.time{trial}(1:numTimePointsToAnalyze);
    else
        error('Trial %d has less than %d time points.', trial, numTimePointsToAnalyze);
    end
end

%% fft that sunuvagun


% alpha
% Prepare the configuration for frequency analysis
cfg = [];
cfg.method = 'mtmfft';       % Multitaper Fourier Transform
cfg.output = 'pow';          % Output the power spectrum
cfg.taper = 'hanning';       % Using Hanning taper
cfg.foi = 10:1:30;           % Frequencies 
cfg.pad = 'maxperlen';       % Padding, if necessary
cfg.tapsmofrq = 2;  % Smoothing of 2 Hz after Mendoza-Halliday et al.

% Run the frequency analysis
freqAnalysisResult = ft_freqanalysis(cfg, data);

% The output, freqAnalysisResult, will have the power spectra
% for the frequencies of interest across all trials and channels
% Access power as: freqAnalysisResult.powspctrm

% Example to display some results
disp('Frequency analysis completed.');
disp(['Number of frequencies analyzed: ' num2str(length(freqAnalysisResult.freq))]);
disp(['Power spectral data size: ' mat2str(size(freqAnalysisResult.powspctrm))]);


% Assuming freqAnalysisResult.powSpectrum holds the power spectrum data
% with dimensions 21 (channels) x 7 (frequencies)

% Initialize the normalized power spectrum array
normalizedPowSpectrum = zeros(size(freqAnalysisResult.powspctrm));

% Loop over each frequency
for freqIdx = 1:size(freqAnalysisResult.powspctrm, 2)
    % Extract the power values for all channels at the current frequency
    powerAtFreq = freqAnalysisResult.powspctrm(:, freqIdx);
    
    % Find the maximum power at this frequency
    maxPower = max(powerAtFreq);
    
    % Normalize the power values by the maximum power at this frequency
    normalizedPowSpectrum(:, freqIdx) = powerAtFreq / maxPower;
end

% The variable 'normalizedPowSpectrum' now contains the normalized power values
% where each power value is divided by the maximum power at that frequency across all channels

% Display or further process 'normalizedPowSpectrum' as needed
disp('Normalized Power Spectrum:');
disp(normalizedPowSpectrum);

averagePowerPerFreq = mean(normalizedPowSpectrum, 2); 
plot(averagePowerPerFreq,1:numChannels)
set(gca, 'YDir', 'reverse');

%% gamma
% Prepare the configuration for frequency analysis
cfg = [];
cfg.method = 'mtmfft';       % Multitaper Fourier Transform
cfg.output = 'pow';          % Output the power spectrum
cfg.taper = 'hanning';       % Using Hanning taper
cfg.foi = 50:150;           % Frequencies
cfg.pad = 'maxperlen';       % Padding, if necessary
cfg.tapsmofrq = 2;  % Smoothing of 2 Hz

% Run the frequency analysis
freqAnalysisResult = ft_freqanalysis(cfg, data);

% The output, freqAnalysisResult, will have the power spectra
% for the frequencies of interest across all trials and channels
% Access power as: freqAnalysisResult.powspctrm

% Example to display some results
disp('Frequency analysis completed.');
disp(['Number of frequencies analyzed: ' num2str(length(freqAnalysisResult.freq))]);
disp(['Power spectral data size: ' mat2str(size(freqAnalysisResult.powspctrm))]);


% Assuming freqAnalysisResult.powSpectrum holds the power spectrum data
% with dimensions 21 (channels) x 7 (frequencies)

% Initialize the normalized power spectrum array
normalizedPowSpectrum = zeros(size(freqAnalysisResult.powspctrm));

% Loop over each frequency
for freqIdx = 1:size(freqAnalysisResult.powspctrm, 2)
    % Extract the power values for all channels at the current frequency
    powerAtFreq = freqAnalysisResult.powspctrm(:, freqIdx);
    
    % Find the maximum power at this frequency
    maxPower = max(powerAtFreq);
    
    % Normalize the power values by the maximum power at this frequency
    normalizedPowSpectrum(:, freqIdx) = powerAtFreq / maxPower;
end

% The variable 'normalizedPowSpectrum' now contains the normalized power values
% where each power value is divided by the maximum power at that frequency across all channels

% Display or further process 'normalizedPowSpectrum' as needed
disp('Normalized Power Spectrum:');
disp(normalizedPowSpectrum);

averagePowerPerFreq = mean(normalizedPowSpectrum, 2); 
hold on
plot(averagePowerPerFreq,1:numChannels)
set(gca, 'YDir', 'reverse');
legend('A-B','Broad Gamma')


