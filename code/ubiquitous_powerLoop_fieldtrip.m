% find spectrolaminar trends using field trip
% chase m 2024

%% Initialization
tic;
clear;
clearvars;
close all;

%% Set directory and parameters
mydir = 'E:\spectrolaminar\AttnData\belt\cont\'; % Specify directory
%mydir = 'E:\spectrolaminar\AttnData\core\cont\'; % Specify directory
figuresDir = fullfile(mydir, 'Crosses_taper'); % dir where the figs go
myfiles = dir(fullfile(mydir,'*@os*.mat')); % Get all files in struct
%myfiles = dir(fullfile(mydir,'*@oe*.mat')); % Get all files in struct
peterdataflag = 1; 
% Frequency bands to analyze
frequencyBands = {
    'Alpha-Beta', 10:30;
    'Broadband Gamma', 50:150;
};

for loopct = 1:length(myfiles)
    %% data3D with dimensions channels, trials, epoched time
    % Load the cont file
    basefilename = myfiles(loopct).name;
    fullfilename = fullfile(mydir, basefilename);
    
    %% is this a peter file or an annie file?
    if peterdataflag ==1
        [~, ~, ~, ~, ~,~,~,lfp] = EphysExtractFxn(fullfilename,1);
    else
        load(fullfilename);
        newadrate = 1000;
        epoch_tframe = [-300 300]; 
        trig0 = anatrig{1};
        trig01 = [];
        for trigredxct=1:length(trig0)
            trig01(trigredxct)    = round(trig0(trigredxct)./(adrate/newadrate));
        end
        
        [lfp] = EphysEpochFxnSimple(cnt,trig01,epoch_tframe,newadrate);

    end
    %% reorg for field trip

    % LFP is a 3D chan x trial x time array sampled at 1 khz
    data3D = lfp;
    data = struct();
    data.label = {};      % Update this with your channel names if needed
    data.fsample = 1000;  % sampling frequency

    % Number of channels, trials, and time points
    [numChannels, numTrials, numTimePoints] = size(data3D);

    % Set channel labels (assuming generic labels here)
    data.label = arrayfun(@(x) sprintf('Ch%d', x), 1:numChannels, 'UniformOutput', false);

    % Reshape 3D double array data for field trip
    data.trial = cell(1, numTrials);
    data.time = cell(1, numTrials);

    % Populate the trial and time fields
    for trial = 1:numTrials
            data.trial{trial} = squeeze(data3D(:, trial, :));
            data.time{trial} = (0:numTimePoints-1) / data.fsample;
    end

    % Want to only analyze part of the epoch?
    numTimePointsToAnalyze = 200; % Number of time points to analyze
    for trial = 1:length(data.trial)
        if length(data.time{trial}) >= numTimePointsToAnalyze
            data.trial{trial} = data.trial{trial}(:, 1:numTimePointsToAnalyze);
            data.time{trial} = data.time{trial}(1:numTimePointsToAnalyze);
        else
            error('Trial %d has less than %d time points.', trial, numTimePointsToAnalyze);
        end
    end

    % Initialize figure for plotting
    figure;
    hold on;

    % Loop through each frequency band
    for bandIdx = 1:size(frequencyBands, 1)
        bandName = frequencyBands{bandIdx, 1};
        freqRange = frequencyBands{bandIdx, 2};

        % Prepare the configuration for frequency analysis
        cfg = [];
        cfg.method = 'mtmfft';       % Multitaper Fourier Transform
        cfg.output = 'pow';          % Output the power spectrum
        cfg.taper = 'hanning';       % Using Hanning taper
        cfg.foi = freqRange;         % Frequencies of interest
        cfg.pad = 'maxperlen';       % Padding, if necessary
        cfg.tapsmofrq = 2;           % Smoothing of 2 Hz

        % Run the frequency analysis
        freqAnalysisResult = ft_freqanalysis(cfg, data);

        % Normalize the power spectrum
        normalizedPowSpectrum = zeros(size(freqAnalysisResult.powspctrm));
        for freqIdx = 1:size(freqAnalysisResult.powspctrm, 2)
            powerAtFreq = freqAnalysisResult.powspctrm(:, freqIdx);
            maxPower = max(powerAtFreq);
            normalizedPowSpectrum(:, freqIdx) = powerAtFreq / maxPower;
        end

        % Calculate average power per frequency
        averagePowerPerFreq = mean(normalizedPowSpectrum, 2);

        % Plot the results
        plot(averagePowerPerFreq, 1:numChannels, 'LineWidth', 2, 'DisplayName', bandName);
    end

    set(gca, 'YDir', 'reverse');
    legend;
    xlabel('Normalized Power');
    ylabel('Channel');
    title(['Power Spectrum: ', basefilename(1:13)]);

    %% Save the figure
    if ~exist(figuresDir, 'dir')
        mkdir(figuresDir);
    end

    % Save the figure with the correct file name
    figureFileName = fullfile(figuresDir, [basefilename(1:end-4), '_PowerSpectrum.fig']);
    saveas(gcf, figureFileName);
    % Save as .jpg
    saveas(gcf, [figureFileName, '.jpg']);

    close(gcf);

    disp(['Processed and saved: ', basefilename]);
end

toc;
