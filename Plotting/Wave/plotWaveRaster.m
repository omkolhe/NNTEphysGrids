function [] = plotWaveRaster(Waves1,Waves2,behaviourTrace1,behaviourTrace2,parameters)

if ~isempty(Waves2)
    waves2Flag = 1;
else
    waves2Flag = 0;
end

waves1Present = vertcat(Waves1.waveStart);
if waves2Flag == 1
    waves2Present = vertcat(Waves2.waveStart);
end

figure();
if waves2Flag==1 ax1=subplot(4,1,1), else ax1=subplot(2,1,1), end
rasterPlot(waves1Present);hold on;
% xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
% plot(RTTraceTime,1:size(IntanBehaviour1.cueHitTrace,2),'.r');
% xline((mean(IntanBehaviour1.reactionTime,'all')*parameters.Fs + parameters.windowBeforeCue*parameters.Fs+1),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([0 size(waves1Present,2)]);%ylim([1 20]);
set(gca,'TickDir','out','fontsize',14'); box off;

if waves2Flag==1 ax2=subplot(4,1,2), else ax2=subplot(2,1,2), end
bar((sum(waves1Present,1)/size(behaviourTrace1,2)));
% xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
% xline((mean(IntanBehaviour1.reactionTime,'all')*parameters.Fs + parameters.windowBeforeCue*parameters.Fs+1),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylim([0 0.4]);xlim([0 size(waves1Present,2)]);
ylabel('Wave probability');xlabel('Time (in ms)');set(gca,'TickDir','out','fontsize',14'); box off;
if waves2Flag == 1
    ax3 = subplot(4,1,3);
    rasterPlot(waves2Present);hold on;
%     xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
    ylabel('Trials');xlabel('Time (in ms)');xlim([0 size(waves2Present,2)]);%ylim([1 20]);
    set(gca,'TickDir','out','fontsize',14'); box off;
    ax4 = subplot(4,1,4)
    bar((sum(waves2Present,1)/size(behaviourTrace2,2)));
    xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
    ylim([0 0.4]);xlim([0 size(waves2Present,2)])
    ylabel('Wave probability');xlabel('Time (in ms)')
    set(gca,'TickDir','out','fontsize',14'); box off;
end

if waves2Flag == 1
    linkaxes([ax1 ax2 ax3 ax4],'x');
else
    linkaxes([ax1 ax2],'x');
end
