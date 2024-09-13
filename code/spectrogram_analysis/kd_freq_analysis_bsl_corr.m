%% Spectrogram ob 30-150 Hz band

% inputs 
% - data: lfp or csd, as cell array of length trials
% - winl: window length in CYCLES
% - correct_1f: boolean: perform 1/f correction?
% - pre_stim: pre-stimulus interval
% - stim: stimulus interval
% - post_stim: post-stimulus interval
% 
% 
% output
% - tfr_lowfreq: spectrogram


function [new_tfr,tfr_pre, tfr_stim, tfr_post] = kd_freq_analysis_bsl_corr(data,stgs)

    % spectrogram
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.pad = 'nextpow2';                   % zero pad to next highest power of 2
    cfg.foi = stgs.freqvec;
    cfg.output = 'pow';
    cfg.toi = -1.75:0.05:1.75;            % center over these time points
    cfg.keeptrials = stgs.keeptrials;

    if stgs.freqvec(1) >= 30
        cfg.t_ftimwin  = stgs.winl./cfg.foi;
    elseif stgs.freqvec(1) < 30
        cfg.t_ftimwin = ones(size(cfg.foi)).*stgs.winl;
    end
    
    tfr = ft_freqanalysis(cfg,data);
    
    %% 1/f correct by means of baseline correction
    
    if stgs.correct_1f

        tiledlayout(2,1)
        nexttile
        plot(tfr.freq,mean(squeeze(mean(tfr.powspctrm,1)),3))
        title('pre bsl corr')
        cfg = [];
        cfg.baselinetype = 'relchange';
        % for low frequencies, use -1 to 1 sec interval (correction to
        % mean)
        if stgs.freqvec(1) < 30
            cfg.baseline = [-1 1];

        % for gamma band + use pre-stimulus interval
        elseif stgs.freqvec(1) >= 30
            cfg.baseline = [-0.75 -.25];
        end
        new_tfr = ft_freqbaseline(cfg,tfr);

        nexttile
        plot(new_tfr.freq,mean(squeeze(mean(new_tfr.powspctrm,1)),3))
        title('post bsl corr')
    end
    
    % split into pre, stimulus onset, and post
    cfg = [];
    cfg.latency = stgs.pre_stim;
    tfr_pre = ft_selectdata(cfg,new_tfr);
    cfg.latency = stgs.stim;
    tfr_stim = ft_selectdata(cfg,new_tfr);
    cfg.latency = stgs.post_stim;
    tfr_post = ft_selectdata(cfg,new_tfr);