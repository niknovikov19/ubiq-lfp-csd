% Input and output folders
dirpath_in = 'D:\WORK\Salvador\repo\ubiquitous-main';
dirpath_out = 'D:\WORK\Salvador\repo\ubiquitous-main\lfp_csd\results';

% % Load LFP data
% data_type = 'LFP';
% X = load(fullfile(dirpath_in, 'RelPowAreaA1.mat'));
% X = X.RelPow;

% Load CSD data
data_type = 'CSD';
X = load(fullfile(dirpath_in, 'RelPowAreaA1_CSD.mat'));
X = X.results;

% Load info about L4 sink channels
Y = load('RelPowAreaA1.mat');
l4_chans = [Y.RelPow.L4chan];
if strcmp(data_type, 'CSD')
    l4_chans(5) = [];  % this lfp recording is missing in csd data
end

for n = 1 : length(X)
    fprintf(num2str(n));
    % Plot
    figure(111); clf; hold on;
    W = X(n).relpow;
    imagesc(W);
    caxis([0, 1]);
    plot([1, nfreq], l4_chans(n) * [1, 1], 'k--');
    xlim([1, size(W, 2)]);
    ylim([1, size(W, 1)]);
    xlabel('Frequency');
    ylabel('Channel');
    set(gca, 'YDir', 'reverse');
    title(sprintf('%s (%s)', X(n).filename, data_type));
    colorbar;
    % Save
    fname_out = sprintf('%i.png', n);
    fpath_out = fullfile(dirpath_out, data_type, fname_out);
    saveas(gcf, fpath_out)
end
