% find spectrolaminar trends using field trip and FLIP analysis
% chase m 2024

%% Initialization
tic;
clear;
clearvars;
close all;

%% Set directory and parameters
mydir = 'E:\spectrolaminar\AttnData\core\cont_offBF\'; % Specify directory
figuresDir = fullfile(mydir, 'Crosses_taper_vFLIPv3'); % dir where the figs go
allFlipResults = struct();
pickchannels = 0; % do you want to manually select channels based on the CSD? 
loadchannels = 0; % do you want to load previously selected channels?
peterdataflag = 1; % peter and noelle's A1 and belt Data

if loadchannels ==1
    % dir with flip results and channel selections
    SelChansDir = fullfile(mydir,'Crosses_taper_FLIP'); 
end

% handle file extension differences
if peterdataflag == 1
    myfiles = dir(fullfile(mydir, '*@os*.mat')); % Get all aud files in struct
else
    myfiles = dir(fullfile(mydir, '*@oe*.mat')); % annie's vis data
end

if size(myfiles, 1) == 0
    disp('no files found!')
end

% Load AllFlipResults if it exists because we want to use our saved
% "Selected Channels variable" for the visual cortex data only!
if loadchannels ==1
    resultsFile = fullfile(SelChansDir, 'allFlipResults.mat');
    if exist(resultsFile, 'file')
        load(resultsFile, 'allFlipResults');
    else
        allFlipResults = struct();
    end
end

for loopct = 1:1%length(myfiles)
    close all
    %% data3D with dimensions channels, trials, epoched time
    % Load the cont file
    basefilename = myfiles(loopct).name;
    fullfilename = fullfile(mydir, basefilename);
    
    % olive recordings had 200 micron spacing
    if peterdataflag == 0
        if basefilename(3) == 'O'
            interelectrodespacing = 0.2;
        else
            interelectrodespacing = 0.1;
        end
    else
        interelectrodespacing = 0.1;
    end
    
    %% is this a peter file or an annie file?
    if peterdataflag == 1
        [~, ~, ~, ~, ~,~,~,lfp] = EphysExtractFxn(fullfilename, 1);
        
    else
        load(fullfilename);
        newadrate = 1000;
        epoch_tframe = [-300 300]; 
        trig0 = anatrig{1};
        trig01 = [];
        for trigredxct = 1:length(trig0)
            trig01(trigredxct) = round(trig0(trigredxct) / (adrate / newadrate));
        end
        [lfp] = EphysEpochFxnSimple(cnt, trig01, epoch_tframe, newadrate);
        
    end
    
    %% Check if selected channels exist in allFlipResults
    if loadchannels ==1
        if isfield(allFlipResults, 'filename') && any(strcmp({allFlipResults.filename}, basefilename))
            currentFileIndex = find(strcmp({allFlipResults.filename}, basefilename));
            if isfield(allFlipResults(currentFileIndex).results, 'selectedChannels')
                selectedChannels = allFlipResults(currentFileIndex).results.selectedChannels;
                lfp = lfp(selectedChannels,:,:);
            end
        end
    end

    %% Select channels (optional)
    if pickchannels == 1
        % Number of channels, trials, and time points
        [numChannels, ~, ~] = size(lfp);

        % Ask user to select channels
        fprintf('Processing file: %s\n', basefilename);
        prompt = sprintf('Select channels (1-%d) separated by commas: ', numChannels);
        selectedChannels = input(prompt, 's');
        selectedChannels = str2num(selectedChannels);  % Convert string to numeric array

        % Validate selected channels
        if isempty(selectedChannels) || any(selectedChannels < 1) || any(selectedChannels > numChannels)
            error('Invalid channel selection.');
        end

        % Extract selected channels
        lfp = lfp(selectedChannels, :, :);
        numChannels = length(selectedChannels); % Update number of channels
    end
    
    %% reorg for field trip
    % LFP is a 3D chan x trial x time array sampled at 1 kHz
    data3D = lfp;
    data = struct();
    data.fsample = 1000; % sampling frequency

    % Number of channels, trials, and time points
    [numChannels, numTrials, numTimePoints] = size(data3D);

    
    % input args for the Mendoza-Halliday, Major et al. code
    freqaxis = 1:150; % put in frequencies of interest
    setfreqbool = 0;
    %data3Dreshaped = reshape(data3D, [numChannels, numTimePoints, numTrials]);
    data3Dreshaped = [];
    for trial_index = 1:numTrials
      data3Dreshaped(:,:,trial_index) = lfp(:,trial_index,:);
    end
    
    % relative power and avg. power calcualation
    [relpow1, nonnormpow1] = relpow_from_rawLFP(data3Dreshaped);
    freqaxis = 1:size(relpow1, 2); % needs to match the freqs from the fft
    laminaraxis = 0:interelectrodespacing:numChannels*interelectrodespacing-interelectrodespacing;
    probe1 = nonnormpow1;

    % run FLIP analysis
    if size(probe1, 1) == 0
        goodnessvalue = NaN;
        lowfreqm = NaN;
        highfreqm = NaN;
        superficialchan = NaN;
        deepchan = NaN;
    else
        [startinglowfreq, endinglowfreq, startinghighfreq, endinghighfreq, goodnessvalue, superficialchannel, deepchannel, highfreqmaxchannel, lowfreqmaxchannel, crossoverchannel] = ...
            FLIPAnalysis(probe1, laminaraxis, freqaxis, setfreqbool);
        
        FlipResults = struct();
        FlipResults.G = goodnessvalue;
        FlipResults.supra = highfreqmaxchannel;
        FlipResults.infra = lowfreqmaxchannel;
        FlipResults.cross = crossoverchannel;
        FlipResults.startinglowfreq = startinglowfreq;
        FlipResults.endinglowfreq = endinglowfreq;
        FlipResults.endinglowfreq = endinglowfreq;
        FlipResults.startinghighfreq = startinghighfreq;
        FlipResults.endinghighfreq = endinghighfreq;
        
        if pickchannels == 1
            FlipResults.selectedChannels = selectedChannels; % Save selected channels
        end
        
        
        allFlipResults(loopct).filename = basefilename;
        allFlipResults(loopct).results = FlipResults;

        if isempty(lowfreqmaxchannel)
            lowfreqm = NaN;
        else
            lowfreqm = lowfreqmaxchannel;
        end
        if isempty(highfreqmaxchannel)
            highfreqm = NaN;
        else
            highfreqm = highfreqmaxchannel;
        end
        if isempty(superficialchannel)
            superficialchan = NaN;
        else
            superficialchan = superficialchannel;
        end
        if isempty(deepchannel)
            deepchan = NaN;
        else
            deepchan = deepchannel;
        end
    end

    %% Save the figure
    if ~exist(figuresDir, 'dir')
        mkdir(figuresDir);
    end

    % Check if a figure was generated by the FLIPAnalysis function
    figHandles = findall(0, 'Type', 'figure');
    if ~isempty(figHandles)
        % Save the figure with the correct file name
        figureFileName = fullfile(figuresDir, [basefilename(1:end-4), '_FLIPAnalysis.fig']);
        saveas(gcf, figureFileName);
        % Save as .jpg
        saveas(gcf, [figureFileName(1:end-4), '.jpg']);
        close(gcf);
    end

    disp(['Processed and saved: ', basefilename]);
end

if ~exist(figuresDir, 'dir')
    mkdir(figuresDir);
end

% Save all FlipResults to a MAT file
save(fullfile(figuresDir, 'allFlipResults.mat'), 'allFlipResults');

toc;
