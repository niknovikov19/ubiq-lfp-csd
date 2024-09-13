%% Make nested plot

%% Inputs
% - chan: channel index
% - data_highfreq, data_lowfreq: struct, output of ft_freqanalysis
% - highfreq_pre, lowfreq_pre: result of ft_freqanalysis for high and low frequencies, pre-stimulus
% - highfreq_stim, lowfreq_stim: same as above, during stimulus
% - highfreq_post, lowfreq_post: same as above post stimulus

function make_nested_plot(chan,data_highfreq,data_lowfreq,stgs,cmap)

    % split into pre, stimulus onset, and post
    cfg = [];
    cfg.latency = stgs.pre_stim;
    cfg.avgovertime = 'yes';
    cfg.avgoverrpt = 'yes';
    highfreq_pre = ft_selectdata(cfg,data_highfreq);
    lowfreq_pre = ft_selectdata(cfg,data_lowfreq);
    
    cfg.latency = stgs.stim;
    highfreq_stim = ft_selectdata(cfg,data_highfreq);
    lowfreq_stim = ft_selectdata(cfg,data_lowfreq);
    
    cfg.latency = stgs.post_stim;
    highfreq_post = ft_selectdata(cfg,data_highfreq);
    lowfreq_post = ft_selectdata(cfg,data_lowfreq);

    % avg over trials
    cfg = [];
    cfg.avgoverrpt = 'yes';
    data_highfreq = ft_selectdata(cfg,data_highfreq);
    data_lowfreq = ft_selectdata(cfg,data_lowfreq);

    t = tiledlayout(2,3,'TileSpacing','compact');
    backgr_ax = axes(t,'XTick',[],'YTick',[],'Box','off');
    backgr_ax.Layout.TileSpan = [2 1];

    %% high
    % spectrogram high frequency
    ax1 = axes(t);
    imagesc(ax1,data_highfreq.time,data_highfreq.freq,squeeze(data_highfreq.powspctrm(chan,:,:)))
    caxis([-2 2])
    ax1.Layout.Tile = 1;
    axis xy
    ax1.Box = 'off';
    ax1.XTick = [];
    cb = colorbar;
    hold on
    plot(zeros(1,1000),linspace(data_highfreq.freq(1),data_highfreq.freq(end),1000),"Color",[0,0,0],'LineStyle','-.','LineWidth',2)


    % spectrum pre high frequency
    ax2 = axes(t);
    bar_plot = bar(ax2,data_highfreq.freq,mean(highfreq_pre.powspctrm(chan,:,:),3));
    bar_plot.Horizontal = 'on';
    ax2.Layout.Tile = 2;
    title('pre')
    ax2.Box = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    
    % spectrum evoked high frequency
    ax3 = axes(t);
    bar_plot = bar(ax3,data_highfreq.freq,mean(highfreq_stim.powspctrm(chan,:,:),3));
    bar_plot.Horizontal = 'on';
    ax3.Layout.Tile = 3;
    title('stim')
    ax3.Box = 'off';
    ax3.XTick = [];
    ax3.YTick = [];

    ax2.XLim = [-1 1];
    ax3.XLim = [-1 1];

    %% low
    % spectrogram 
    ax4 = axes(t);
    imagesc(ax4,data_lowfreq.time,data_lowfreq.freq,squeeze(data_lowfreq.powspctrm(chan,:,:)))
    caxis([-2 2])
    ax4.Layout.Tile = 4;
    axis xy
    ax4.Box = 'off';
    cb = colorbar;
    hold on

    plot(zeros(1,1000),linspace(data_lowfreq.freq(1),data_lowfreq.freq(end),1000),"Color",[0,0,0],'LineStyle','-.','LineWidth',2)


    % spectrum pre 
    ax5 = axes(t);
    bar_plot = bar(ax5,data_lowfreq.freq,mean(lowfreq_pre.powspctrm(chan,:,:),3));
    bar_plot.Horizontal = 'on';
    ax5.Layout.Tile = 5;
    ax5.YTick = [];
    ax5.Box = 'off';
    % spectrum evoked
    ax6 = axes(t);
    bar_plot = bar(ax6,data_lowfreq.freq,mean(lowfreq_stim.powspctrm(chan,:,:),3));
    bar_plot.Horizontal = 'on';
    ax6.Layout.Tile = 6;
    ax6.YTick = [];
    ax6.Box = 'off';
    ax5.XLim = [-1 1];
    ax6.XLim =[-1 1];

    title(t,data_lowfreq.label{chan},'FontSize', 24)

    ax1.FontSize = 20;
    ax2.FontSize = 20;
    ax3.FontSize = 20;
    ax4.FontSize = 20;
    ax5.FontSize = 20;
    ax6.FontSize = 20;

    colormap(cmap)