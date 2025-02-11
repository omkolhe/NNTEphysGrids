%% Plotting CCA for each session 
argIn = combinedCCA(1).argIn;
mapDim = size(combinedCCA(1).CCA.cueHit.CorrMap, 2);
delays = (-argIn.MaxDelay:argIn.MaxDelay); % Convert to ms
t = (argIn.WindowLength/2:argIn.TimeStep:argIn.WindowLength/2+argIn.TimeStep*(mapDim-1)); % Convert to ms
zeroDelayIndex = floor(numel(delays)/2) + 1;

figure();sgtitle('CCA - Hits vs Misses')
for CANONICAL_PAIR_IDX = 1:3
    avgCCA.cueHit = cell2mat(arrayfun(@(s) squeeze(s.CCA.cueHit.CorrMap(zeroDelayIndex,:,CANONICAL_PAIR_IDX))', combinedCCA, 'UniformOutput', false));
    avgCCA.cueMiss = cell2mat(arrayfun(@(s) squeeze(s.CCA.cueMiss.CorrMap(zeroDelayIndex,:,CANONICAL_PAIR_IDX))', combinedCCA, 'UniformOutput', false));
    
    subplot(3,1,CANONICAL_PAIR_IDX)
    h1=plot(t,squeeze(mean(avgCCA.cueHit,2)),'Color', [0 0.1 0.8],'LineWidth',2); hold on;
    plot(t,squeeze(mean(avgCCA.cueHit,2)) - (squeeze(std(avgCCA.cueHit,1,2)))/sqrt(size(avgCCA.cueHit,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
    plot(t,squeeze(mean(avgCCA.cueHit,2)) + (squeeze(std(avgCCA.cueHit,1,2)))/sqrt(size(avgCCA.cueHit,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
    h2=plot(t,squeeze(mean(avgCCA.cueMiss,2)),'Color', [0.5 0.5 0.5],'LineWidth',2);
    plot(t,squeeze(mean(avgCCA.cueMiss,2)) - (squeeze(std(avgCCA.cueMiss,1,2)))/sqrt(size(avgCCA.cueMiss,2)) ,'Color', [0.5 0.5 0.5 0.4],'LineWidth',2);
    plot(t,squeeze(mean(avgCCA.cueMiss,2)) + (squeeze(std(avgCCA.cueMiss,1,2)))/sqrt(size(avgCCA.cueMiss,2)) ,'Color', [0.5 0.5 0.5 0.4],'LineWidth',2);
    xline(1501,'--r','Cue');%xline(1500+parameters.Fs*mean(IntanBehaviour.reactionTime,'all'),'--r','RT');
    xlabel('Time (ms)'); ylabel('Population Correlation');
    legend([h1 h2],'Hit','Miss','Location','best'); %ylim([5 15]);
    title(['Dimension - ',string(CANONICAL_PAIR_IDX)])
    xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');
    drawnow;
end

% figure();sgtitle('CCA - Hits vs Misses')
% for CANONICAL_PAIR_IDX = 1:3
%     avgCCA.cueHit = cell2mat(arrayfun(@(s) squeeze(s.CCA.cueHit.CorrMap(zeroDelayIndex,:,CANONICAL_PAIR_IDX))', combinedCCA, 'UniformOutput', false));
%     
%     subplot(1,3,CANONICAL_PAIR_IDX)
%     h1=plot(t,avgCCA.cueHit,'LineWidth',2); hold on;
%     xline(1501,'--r','Cue');%xline(1500+parameters.Fs*mean(IntanBehaviour.reactionTime,'all'),'--r','RT');
%     xlabel('Time (ms)'); ylabel('Population Correlation');
%     xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');
%     drawnow;
% end

figure(); sgtitle("CCA - MI Hits vs MI FAs")
for CANONICAL_PAIR_IDX = 1:3
    avgCCA.MIHit = cell2mat(arrayfun(@(s) squeeze(s.CCA.MIHit.CorrMap(zeroDelayIndex,:,CANONICAL_PAIR_IDX))', combinedCCA, 'UniformOutput', false));
    avgCCA.MIFA = cell2mat(arrayfun(@(s) squeeze(s.CCA.MIFA.CorrMap(zeroDelayIndex,:,CANONICAL_PAIR_IDX))', combinedCCA, 'UniformOutput', false));
    
    subplot(3,1,CANONICAL_PAIR_IDX)
    h1=plot(t,squeeze(mean(avgCCA.MIHit,2)),'Color', [0 0.1 0.8],'LineWidth',2); hold on;
    plot(t,squeeze(mean(avgCCA.MIHit,2)) - (squeeze(std(avgCCA.MIHit,1,2)))/sqrt(size(avgCCA.MIHit,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
    plot(t,squeeze(mean(avgCCA.MIHit,2)) + (squeeze(std(avgCCA.MIHit,1,2)))/sqrt(size(avgCCA.MIHit,2)) ,'Color', [0 0.1 0.8 0.4],'LineWidth',2);
    h2=plot(t,squeeze(mean(avgCCA.MIFA,2)),'Color', [0.9 0.1 0.1],'LineWidth',2);
    plot(t,squeeze(mean(avgCCA.MIFA,2)) - (squeeze(std(avgCCA.MIFA,1,2)))/sqrt(size(avgCCA.MIFA,2)) ,'Color', [0.9 0.1 0.1 0.4],'LineWidth',2);
    plot(t,squeeze(mean(avgCCA.MIFA,2)) + (squeeze(std(avgCCA.MIFA,1,2)))/sqrt(size(avgCCA.MIFA,2)) ,'Color', [0.9 0.1 0.1 0.4],'LineWidth',2);
    xline(1501,'--r','MI');%xline(1500+parameters.Fs*mean(IntanBehaviour.reactionTime,'all'),'--r','RT');
    xlabel('Time (ms)'); ylabel('Population Correlation');
    legend([h1 h2],'MIHit','MIFA','Location','best'); %ylim([5 15]);
    title(['Dimension - ',string(CANONICAL_PAIR_IDX)])
    xlim([0 3000]);box off;set(gca,'TickDir','out','fontsize',14');
    drawnow;
end

%% Ploting with the delay 

CANONICAL_PAIR_IDX = 1;
