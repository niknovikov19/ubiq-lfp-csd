% ubiquitous analysis of relative power as a function of brain area
% chase m 2024
% Load the data.mat file
clear; clc; close all;
load('data.mat');

% Extract the relevant fields from the data struct
relpow = data.relpow;
meta = data.meta;

% Reshape the relpow matrix
[num_probes, num_channels, num_freqs] = size(relpow);
relpow_flat = reshape(relpow, [num_probes * num_channels * num_freqs, 1]);

% Create vectors for frequencies and channels
frequencies = repmat(1:num_freqs, [num_probes * num_channels, 1]);
frequencies_flat = reshape(frequencies, [num_probes * num_channels * num_freqs, 1]);

channels = repmat((1:num_channels)', [num_probes, num_freqs]);
channels_flat = reshape(channels, [num_probes * num_channels * num_freqs, 1]);

% Ensure meta.monkey_number is a column vector
monkey_number = [meta.monkey_number]';

% Repeat brain_area, filename, and monkey_number for all channels and frequencies
brain_area_repeated = repmat({meta.brain_area}', [1, num_channels * num_freqs]);
brain_area_flat = reshape(brain_area_repeated, [num_probes * num_channels * num_freqs, 1]);

filename_repeated = repmat({meta.filename}', [1, num_channels * num_freqs]);
filename_flat = reshape(filename_repeated, [num_probes * num_channels * num_freqs, 1]);

monkey_number_repeated = repmat(monkey_number, [1, num_channels * num_freqs]);
monkey_number_flat = reshape(monkey_number_repeated, [num_probes * num_channels * num_freqs, 1]);

% Remove rows with NaN values in relpow
valid_idx = ~isnan(relpow_flat);
relpow_valid = relpow_flat(valid_idx);
frequency_valid = frequencies_flat(valid_idx);
channel_valid = channels_flat(valid_idx);
brain_area_valid = brain_area_flat(valid_idx);
filename_valid = filename_flat(valid_idx);
monkey_number_valid = monkey_number_flat(valid_idx);

% Create a table for the analysis
tbl = table(relpow_valid, frequency_valid, channel_valid, ...
            categorical(brain_area_valid), categorical(filename_valid), ...
            monkey_number_valid, ...
            'VariableNames', {'RelPow', 'Frequency', 'Channel', 'BrainArea', 'Filename', 'MonkeyNumber'});

% Define the formula for the mixed-effects model including BrainArea
formula = 'RelPow ~ Frequency * Channel * BrainArea + (1|Filename) + (1|MonkeyNumber)';

% Fit the mixed-effects model
% Model fit statistics:
%     AIC            BIC           LogLikelihood    Deviance   
%     -6.7779e+05    -6.771e+05    3.3895e+05       -6.7789e+05
lme = fitlme(tbl, formula);

% Display the model summary
disp(lme);

% Perform ANOVA on the fitted model to test for significant effects
anova_tbl = anova(lme);
disp('ANOVA table for the mixed-effects model:');
disp(anova_tbl);


% Calculate residuals
residuals = residuals(lme);

% Plot residuals vs. fitted values
figure;
scatter(predict(lme), residuals);
xlabel('Fitted Values');
ylabel('Residuals');
title('Residuals vs. Fitted Values');

% Plot histogram of residuals
figure;
histogram(residuals);
xlabel('Residuals');
ylabel('Frequency');
title('Histogram of Residuals');



%% modeling brain area as random effect
% Define the formula for the mixed-effects model including BrainArea
formula_AreaRand = 'RelPow ~ Frequency * Channel + (1|BrainArea) + (1|MonkeyNumber)';

% Fit the mixed-effects model
lme_Area_Rand = fitlme(tbl, formula_AreaRand);

% Display the model summary Model fit statistics:
%     AIC      BIC      LogLikelihood    Deviance
%     56790    56885    -28388           56776   
disp(lme_Area_Rand);

% Perform ANOVA on the fitted model to test for significant effects
anova_tbl = anova(lme_Area_Rand);
disp('ANOVA table for the mixed-effects model:');
disp(anova_tbl);
