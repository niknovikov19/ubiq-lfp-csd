%% TFR of LFP, CSD, CSD NNN
clear all; close all; clc

winl_low = 0.5; % window length for low frequenices, in seconds
winl_high = 6; % window length for low frequenices, in cycles

pre_post = [-2, 3.5];           % length of time interval pre- and post, in seconds

fs = 1000;                          % sampling rate

% time intervals for pre- and post-stimulus, and evoked response to
% stimulus (used to get spectra from spectrogram)
pre_stim = [-1. -0.25-1/fs];
stim = [0 0.5-1/fs];
post_stim = [0.5 1.5];

%% Settings
addpath('/Users/katharinaduecker/Documents/MATLAB/fieldtrip')
ft_defaults;

outpth = fullfile('data_exploration');
figpth = fullfile(outpth, 'figures','bslcorr');
mkdir(outpth)
mkdir(figpth)
load(fullfile('data.mat'))

addpath('cbrewer2')
cmap = cbrewer2('div','RdBu',1001);
cmap = flip(cmap);

%% structure LFP in fieldtrip format
LFP = squeeze(data.example1_vlPFC_lfp);
LFP_trials = cell(1,size(LFP,3));
for t = 1:size(LFP,3)
    LFP_trials{t} = LFP(:,:,t);
end


CSD = calc_csd(LFP,1,1);    % CSD

spat_fs = 6666.67; % spatial sampling frequency for 150 micron spaacing

CSDnnn = calc_csd(LFP,2,spat_fs); % CSD, next nearest neighbour

% convert to structure

lfp_trl = cell_to_struct(LFP_trials,1,1000,pre_post);
csd_trl = cell_to_struct(CSD,1,1000,pre_post);
csdnnn_trl = cell_to_struct(CSDnnn,1,1000,pre_post);

%% LFP

% define settings struct (akin to cfg struct in fieldtrip
stgs = [];
stgs.winl = winl_low;
stgs.correct_1f = true;
stgs.pre_stim = pre_stim;
stgs.stim = stim;
stgs.post_stim = post_stim;
stgs.keeptrials = 'yes';                     % keep trials? (no: average)
% 
% % Spectrograms low frequency
% stgs.freqvec = 4:1/winl_low:30;
% fig = figure;
% [lfp_tfr_low, lfp_tfr_low_pre, lfp_tfr_low_stim, lfp_tfr_low_post] = kd_freq_analysis_bsl_corr(lfp_trl,stgs);
% print(fig,fullfile(figpth,'bslcorr_results_lowfreq_all'),'-dpng')
% 
% Spectrograms high frequency
stgs.freqvec = 30:2:150;
stgs.winl = winl_high;
% 
% fig = figure; 
% [lfp_tfr_high, lfp_tfr_high_pre, lfp_tfr_high_stim, lfp_tfr_high_post] = kd_freq_analysis_bsl_corr(lfp_trl,stgs);
% print(fig,fullfile(figpth,'bslcorr_results_highfreq_all'),'-dpng')

% save(fullfile(outpth,'lowfreq_results_manual1f.mat'),'lfp_tfr_low', 'lfp_tfr_low_pre', 'lfp_tfr_low_stim', 'lfp_tfr_low_post')
% save(fullfile(outpth,'highfreq_results_manual1f.mat'),'lfp_tfr_high', 'lfp_tfr_high_pre', 'lfp_tfr_high_stim', 'lfp_tfr_high_post')

load(fullfile(outpth,'lowfreq_results_manual1f.mat'))
load(fullfile(outpth,'highfreq_results_manual1f.mat'))

% plot
for h = 1:length(lfp_tfr_high.label)
    fig = figure('Position',[0 0 1500 500]);
    make_nested_plot(h,lfp_tfr_high,lfp_tfr_low,stgs,cmap)
    print(fig,fullfile(figpth,['lfp_motif_bslcorr_',lfp_tfr_high.label{h}]),'-dpng')
    print(fig,fullfile(figpth,['lfp_motif_bslcorr_',lfp_tfr_high.label{h}]),'-dsvg')
    close all
end

%% CSD
% Spectrograms low frequency
% stgs.freqvec = 4:1/winl_low:30;
% stgs.winl = winl_low;
% 
% fig = figure;
% [csd_tfr_low, csd_tfr_low_pre, csd_tfr_low_stim, csd_tfr_low_post] = kd_freq_analysis_bsl_corr(csd_trl,stgs);
% print(fig,fullfile(figpth,'csd_bslcorr_results_lowfreq_all'),'-dpng')
% 
% % Spectrograms high frequency
% stgs.freqvec = 30:2:150;
% stgs.winl = winl_high;
% fig = figure; 
% [csd_tfr_high, csd_tfr_high_pre, csd_tfr_high_stim, csd_tfr_high_post] = kd_freq_analysis_bsl_corr(csd_trl,stgs);
% print(fig,fullfile(figpth,'csd_bslcorr_results_highfreq_all'),'-dpng')
% 
% save(fullfile(outpth,'lowfreq_results_manual1f_csd.mat'),'csd_tfr_low', 'csd_tfr_low_pre', 'csd_tfr_low_stim', 'csd_tfr_low_post')
% save(fullfile(outpth,'highfreq_results_manual1f_csd.mat'),'csd_tfr_high', 'csd_tfr_high_pre', 'csd_tfr_high_stim', 'csd_tfr_high_post')

load(fullfile(outpth,'lowfreq_results_manual1f_csd.mat'))
load(fullfile(outpth,'highfreq_results_manual1f_csd.mat'))

% plot
for h = 1:length(csd_tfr_high.label)
    fig = figure('Position',[0 0 1500 500]);
    make_nested_plot(h,csd_tfr_high,csd_tfr_low,stgs,cmap)
    print(fig,fullfile(figpth,['csd_motif_bslcorr_',csd_tfr_high.label{h}]),'-dpng')
    print(fig,fullfile(figpth,['csd_motif_bslcorr_',csd_tfr_high.label{h}]),'-dsvg')
    close all
end


%% CSD NNN
% Spectrograms low frequency
% stgs.freqvec = 4:1/winl_low:30;
% stgs.winl = winl_low;
% fig = figure;
% [csdnnn_tfr_low, csdnnn_tfr_low_pre, csdnnn_tfr_low_stim, csdnnn_tfr_low_post] = kd_freq_analysis_bsl_corr(csd_trl,stgs);
% 
% print(fig,fullfile(figpth,'csd_bslcorr_results_lowfreq_all'),'-dpng')
% 
% % Spectrograms high frequency
% stgs.freqvec = 30:2:150;
% stgs.winl = winl_high;
% fig = figure; 
% [csdnnn_tfr_high, csdnnn_tfr_high_pre, csdnnn_tfr_high_stim, csdnnn_tfr_high_post] = kd_freq_analysis_bsl_corr(csd_trl,stgs);
% print(fig,fullfile(figpth,'csd_bslcorr_results_highfreq_all'),'-dpng')
% 
% save(fullfile(outpth,'lowfreq_results_manual1f_csdnnn.mat'),'csdnnn_tfr_low', 'csdnnn_tfr_low_pre', 'csdnnn_tfr_low_stim', 'csdnnn_tfr_low_post')
% save(fullfile(outpth,'highfreq_results_manual1f_csdnnn.mat'),'csdnnn_tfr_high', 'csdnnn_tfr_high_pre', 'csdnnn_tfr_high_stim', 'csdnnn_tfr_high_post')

load(fullfile(outpth,'lowfreq_results_manual1f_csdnnn.mat'))
load(fullfile(outpth,'highfreq_results_manual1f_csdnnn.mat'))

% plot
for h = 1:length(csdnnn_tfr_high.label)
    fig = figure('Position',[0 0 1500 500]);
    make_nested_plot(h,csdnnn_tfr_high,csdnnn_tfr_low,stgs,cmap)
    print(fig,fullfile(figpth,['csdnnn_motif_bslcorr_',csdnnn_tfr_high.label{h}]),'-dpng')
    print(fig,fullfile(figpth,['csdnnn_motif_bslcorr_',csdnnn_tfr_high.label{h}]),'-dsvg')
    close all
end


%% Correlate alpha power and broadband response


% frequency bands to correlate
low_freq = [8 15];
high_freq = [50 100];
fig = figure('Position',[0 0 1000 200]);

% LFP
cfg = [];
cfg.latency = stim;
cfg.frequency = high_freq;
cfg.avgovertime = 'yes';
cfg.avgoverfreq = 'yes';
avg_lfp_pow_high = ft_selectdata(cfg,lfp_tfr_high_stim);

cfg.latency = pre_stim;
cfg.frequency = low_freq;
avg_lfp_pow_alpha = ft_selectdata(cfg,lfp_tfr_low_pre);

% correlation
[corr_alpha_bbg,p] = corr(avg_lfp_pow_alpha.powspctrm,avg_lfp_pow_high.powspctrm,'type','Spearman');

% BH correction
fdr = 0.25; % false discovery rate
m = numel(p); % number of tests
[~,idx] = sort(p);
p_bh = (idx/m).*fdr;
p_sign = p < p_bh;
p_sign(p_bh==0) = 0.25;
redmap = cbrewer2('seq','Reds');
ax(1) = subplot(131);
imagesc(corr_alpha_bbg, 'AlphaData',p_sign)
colormap(ax(1), cmap)
caxis([-.2,.2])
cb = colorbar;
cb.Ticks = -0.2:0.2:0.2;
cb.Label.String = 'Spearman corr';
axis xy
xlabel('channel index alpha power')
ylabel('channel index bb power')
title('LFP')

% CSD
cfg = [];
cfg.latency = stim;
cfg.frequency = high_freq;
cfg.avgovertime = 'yes';
cfg.avgoverfreq = 'yes';
avg_csd_pow_high = ft_selectdata(cfg,csd_tfr_high_stim);

cfg.latency = pre_stim;
cfg.frequency = low_freq;
avg_csd_pow_alpha = ft_selectdata(cfg,csd_tfr_low_pre);

[corr_alpha_bbg,p] = corr(avg_csd_pow_alpha.powspctrm,avg_csd_pow_high.powspctrm,'type','Spearman');
% BH correction
m = numel(p); % number of tests
[~,idx] = sort(p);
p_bh = (idx/m).*fdr;
p_sign = p < p_bh;
p_sign(p_bh==0) = 0.25;
ax(1) = subplot(132);
imagesc(corr_alpha_bbg, 'AlphaData',p_sign)
colormap(ax(1), cmap)
caxis([-.2,.2])
cb = colorbar;
cb.Ticks = -0.2:0.2:0.2;
cb.Label.String = 'Spearman corr';
axis xy
xlabel('channel index alpha power')
ylabel('channel index bb power')
title('CSD')

% CSD NNN
cfg = [];
cfg.latency = stim;
cfg.frequency = high_freq;
cfg.avgovertime = 'yes';
cfg.avgoverfreq = 'yes';
avg_csdnnn_pow_high = ft_selectdata(cfg,csdnnn_tfr_high_stim);

cfg.latency = pre_stim;
cfg.frequency = low_freq;
avg_csdnnn_pow_alpha = ft_selectdata(cfg,csdnnn_tfr_low_pre);

[corr_alpha_bbg,p] = corr(avg_csdnnn_pow_alpha.powspctrm,avg_csdnnn_pow_high.powspctrm,'type','Spearman');

% BH correction
m = numel(p); % number of tests
[~,idx] = sort(p);
p_bh = (idx/m).*fdr;
p_sign = p < p_bh;
p_sign(p_bh==0) = 0.25;

ax(1) = subplot(133);
imagesc(corr_alpha_bbg,'AlphaData',p_sign)
colormap(ax(1), cmap)
caxis([-.2,.2])
cb = colorbar;
cb.Ticks = -0.2:0.2:0.2;
cb.Label.String = 'Spearman corr';
axis xy
xlabel('channel index alpha power')
ylabel('channel index bb power')
title('CSD NNN')

print(fig,fullfile(figpth,'alpha_bb_corr_all'),'-dsvg')
print(fig,fullfile(figpth,'alpha_bb_corr_all'),'-dpng')

%% CSD time

cfg = [];
cfg.avgoverrpt = 'yes';
CSD_time = ft_selectdata(cfg,csd_trl);

close all
fig = figure;
imagesc(CSD_time.time{1}, 1:14,CSD_time.trial{1})
hold on
plot(zeros(1,1000),linspace(1,14,1000),'Color',[0,0,0],'LineStyle','-.','LineWidth',3)
ylim([1 14])
ylabel('channel index')
xlabel('time (s)')
caxis([-0.02 0.02])
cmap = colormap("jet");
cmap = flip(cmap)
colormap(cmap);
cb = colorbar;
cb.Label.String = 'mV/mm^2';
cb.TickLabels{1} = 'sink';
cb.TickLabels{end} = 'source';

print(fig,fullfile(figpth,'csd_time_l4'),'-dsvg')
print(fig,fullfile(figpth,'csd_time_l4'),'-dpng')