% Input and output folders
dirpath_in = 'D:\WORK\Salvador\repo\ubiquitous-main';
dirpath_out = 'D:\WORK\Salvador\repo\ubiquitous-main\lfp_csd\results\LFP_CSD_3';

mkdir(dirpath_out);

% Load data
X = struct();
X.lfp = load(fullfile(dirpath_in, 'RelPowAreaA1.mat'));
X.lfp = X.lfp.RelPow;
X.csd = load(fullfile(dirpath_in, 'RelPowAreaA1_CSD.mat'));
X.csd = X.csd.results;

% Remove mismatching channels
X.lfp(5) = [];  % this lfp recording is missing in csd data
X.lfp(4) = [];
X.csd(4) = [];

nsess = length(X.lfp);
nchan = size(X.lfp(1).relpow, 1);
nfreq = size(X.lfp(1).relpow, 2);

% L4 sink channels
l4_chans = [X.lfp.L4chan];

data_types = {'lfp', 'csd'};
nx = 3;
ny = 2;

%fbands = struct('low', [10, 30], 'high', [40, 100]);
%fbands = struct('low', [10, 30], 'high', [65, 150]);
fbands = struct('low', [1, 10], 'high', [65, 150]);

% Average over sessions
Wavg = struct();
for m = 1 : 2
    data_type = data_types{m};
    W_ = NaN * ones(nchan, nfreq, nsess);
    for n = 1 : nsess
        W_(:, :, n) = X.(data_type)(n).relpow;
    end
    Wavg.(data_type) = mean(W_, 3);
end

for n = 0 : nsess
    
    fprintf('%i\n', n);
    figure(111); clf; hold on;
    
    for m = 1 : 2
       
        data_type = data_types{m};
        
        if n == 0
            l4_chan = mean(l4_chans);
            W = Wavg.(data_type);
            fname = 'avg';
        else
            Xcur = X.(data_type)(n);
            l4_chan = l4_chans(n);
            W = Xcur.relpow;
            fname = Xcur.filename;
        end
        
        %W(:, 60) = (W(:, 59) + W(:, 61)) / 2;
        
        % Depth-frequency
        subplot(ny, nx, (m-1) * nx + [1, 2]); hold on;
        imagesc(W);
        caxis([0, 1]);
        plot([1, nfreq], l4_chan * [1, 1], 'k--');
        xlim([1, size(W, 2)]);
        ylim([1, size(W, 1)]);
        xlabel('Frequency');
        ylabel('Channel');
        set(gca, 'YDir', 'reverse');
        title(sprintf('%s (%s)', fname, data_type));
        colorbar;
        
        % Depth profile
        subplot(ny, nx, (m-1) * nx + 3); hold on;
        for fb = {'low', 'high'}
        %for fb = {'high'}
            fband = fbands.(fb{:});
            w = mean(W(:, fband(1) : fband(2)), 2);
            plot(w, [1 : length(w)], 'DisplayName', sprintf('%i-%i', fband(1), fband(2)));
        end
        ylim([1, length(w)]);
        set(gca, 'YDir', 'reverse');
        %legend('show');
        title_str = sprintf('%i-%i / %i-%i', fbands.low(1), fbands.low(2),...
            fbands.high(1), fbands.high(2));
        title(title_str);
        %title(sprintf('%i-%i Hz', fband(1), fband(2)))
        
    end
    
    % Save
    if n == 0
        fname_out = 'avg.png';
    else
        fname_out = sprintf('%i.png', n);
    end
    fpath_out = fullfile(dirpath_out, fname_out);
    saveas(gcf, fpath_out)
    
end
