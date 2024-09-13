% Analyze FlipResults
% chase m 2024

%% Load AllFlipResults
figuresDir = 'E:\spectrolaminar\VisCtxData\cont\Crosses_taper_vFLIPv2'; % Directory where the results are saved
resultsFile = fullfile(figuresDir, 'allFlipResults.mat');
load(resultsFile);



%% Extract and Analyze Goodness Values (G)
G_values = [];
for i = 1:length(allFlipResults)
    G_values = [G_values; allFlipResults(i).results.G];
end

% Separate G values into bins and label NaNs as "not sig."
G_isnan = isnan(G_values);
G_values_with_nans = G_values; % Preserve original G values including NaNs
G_values(G_isnan) = -5; % Temporarily set NaNs to -5 for binning

% Define bins for histogram, including one bin for "not sig."
binEdges = [-5, -2:0.1:2]; % Adjust bin edges as needed
binLabels = arrayfun(@(x) sprintf('%.1f', x), binEdges(2:end), 'UniformOutput', false);
binLabels = ['not sig', binLabels];

% Calculate histogram
histCounts = histcounts(G_values, binEdges);

% Only label every 10th tick
numLabels = length(binLabels);
step = 3; 
xticks = 1:step:numLabels;
xticklabels = binLabels(1:step:numLabels);

%% Plot Histogram of G Values
figure;
bar(histCounts);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45); % Angle the x-axis labels by 45 degrees
xlabel('Goodness Value (G)');
ylabel('Count');
title('Histogram of Goodness Values (G)');
grid on;

% Save the histogram figure
histogramFileName = fullfile(figuresDir, 'GoodnessValuesHistogram.fig');
saveas(gcf, histogramFileName);
% Save as .jpg
saveas(gcf, [histogramFileName(1:end-4), '.jpg']);
close(gcf);


% distribution of crossover channels
crossoverChannels = [];
for i = 1:length(allFlipResults)
    crossoverChannels = [crossoverChannels; allFlipResults(i).results.cross];
end

figure;
histogram(crossoverChannels, 'BinEdges', 1:max(crossoverChannels)+1);
xlabel('Crossover Channel');
ylabel('Count');
title('Distribution of Crossover Channels');
grid on;

% Save the crossover channel figure
crossoverFileName = fullfile(figuresDir, 'CrossoverChannelsHistogram.fig');
saveas(gcf, crossoverFileName);
% Save as .jpg
saveas(gcf, [crossoverFileName(1:end-4), '.jpg']);
close(gcf);

%statistics for G values
mean_G = mean(G_values_with_nans(~G_isnan));
std_G = std(G_values_with_nans(~G_isnan));
fprintf('Mean G value: %.2f\n', mean_G);
fprintf('Standard deviation of G values: %.2f\n', std_G);

% Initialize an array to store the differences
differences = [];

% Initialize arrays to store the differences and distances
differences = [];
supra_distances = [];
infra_distances = [];

% Loop through each file in allFlipResults
for i = 1:length(allFlipResults)
    % Determine the electrode distance based on filename
    if allFlipResults(i).filename(3) == 'O'
        elec_dist = 0.2;
    else
        elec_dist = 0.1;
    end
    
    if isfield(allFlipResults(i).results, 'supra') && isfield(allFlipResults(i).results, 'infra')
        supra = allFlipResults(i).results.supra * elec_dist;
        infra = allFlipResults(i).results.infra * elec_dist;
        
        % Store the distances
        if ~isnan(supra)
            supra_distances = [supra_distances; supra];
        end
        if ~isnan(infra)
            infra_distances = [infra_distances; infra];
        end

        % Calculate the absolute difference and store it
        if ~isnan(supra) && ~isnan(infra)
            differences = [differences; abs(supra - infra)];
        end
        
        % Get the filename
        filename = allFlipResults(i).filename;
        disp(['Processing file: ', filename]);
    end
end

% Plot histograms
figure;

% Subplot 1: Histogram of absolute differences
subplot(3,1,3);
histogram(differences, 'BinWidth', 0.1);
xlabel('Absolute Difference between Supra and Infra (mm)');
ylabel('Count');
title('Estimates of Supra/Infra distance');
grid on;

% Subplot 2: Histogram of supra distances
subplot(3,1,1);
histogram(supra_distances, 'BinWidth', 0.1);
xlabel('Supra Distance (mm)');
ylabel('Count');
title('Estimates of Supragranular');
grid on;

% Subplot 3: Histogram of infra distances
subplot(3,1,2);
histogram(infra_distances, 'BinWidth', 0.1);
xlabel('Infra Distance (mm)');
ylabel('Count');
title('Estimates of Infragranular');
grid on;

% Save the histogram figure
histogramFileName = fullfile(figuresDir, 'DifferencesSupraInfraDistancesHistogram.fig');
saveas(gcf, histogramFileName);
% Save as .jpg
saveas(gcf, [histogramFileName(1:end-4), '.jpg']);
close(gcf);


toc;
