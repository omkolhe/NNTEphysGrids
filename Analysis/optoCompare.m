%% Comaparing Behaviour 
[p,h] = ranksum(IntanBehaviourBaseline.reactionTime,IntanBehaviourOpto.reactionTime)
plotBox2(IntanBehaviourBaseline.reactionTime,IntanBehaviourOpto.reactionTime);
% xL=xlim;
% yL=ylim;
% text(0.995*xL(2),0.995*yL(2),['p-val = ' num2str(p)],'HorizontalAlignment','right','VerticalAlignment','top')
ylabel('Reaction Time (s)'); title('M2 -> M1 Opto');subtitle(['p-val = ' num2str(p)]);
xtix = {'Baseline','Opto'}; xtixloc = [1 2]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
set(gca,'TickDir','out','fontsize',14');

%% Comparing PA 
PABaseline = getPA(IntanBehaviourBaseline,0,1,0,parameters,0);
PAOpto = getPA(IntanBehaviourOpto,0,1,0,parameters,0);
% Smooth
figure();
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(squeeze(mean(PABaseline.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[0.7 0.7 0.7],'LineWidth',1.5); hold on;
plot(IntanBehaviourOpto.cueHitTrace(1).time,smooth(squeeze(mean(PAOpto.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
ylabel("Phase Alignment"); xlabel("Time (s)");
xline(0,'--k','Cue','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xlim([-0.5 1.5]);box off;legend('Baseline','Opto');legend('boxoff');set(gca,'TickDir','out','fontsize',14');
title("Phase Alignment : M2 -> M1 Opto");
% Raw
figure();
plot(IntanBehaviourBaseline.cueHitTrace(1).time,squeeze(mean(PABaseline.Hit,[1 2],'omitnan')),'Color',[0.7 0.7 0.7],'LineWidth',1.5); hold on;
plot(IntanBehaviourOpto.cueHitTrace(1).time,squeeze(mean(PAOpto.Hit,[1 2],'omitnan')),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
ylabel("Phase Alignment"); xlabel("Time (s)");
xline(0,'--k','Cue','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xlim([-0.5 1.5]);box off;legend('Baseline','Opto');legend('boxoff');set(gca,'TickDir','out','fontsize',14');
title("Phase Alignment : M2 -> M1 Opto");


%% Comparing PDG 
PGD.avgPGDBaseline = mean(vertcat(WavesBaseline.wavesHit.PGD),1);
PGD.avgPGDOpto = mean(vertcat(WavesOpto.wavesHit.PGD),1);

figure(); hold on;
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(PGD.avgPGDBaseline,50,'sgolay',20),'Color',[0.7 0.7 0.7],'LineWidth',1.5);
plot(IntanBehaviourOpto.cueHitTrace(1).time,smooth(PGD.avgPGDOpto,50,'sgolay',20),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off; legend('Baseline','Opto');legend('boxoff');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');

%% Comparing Wave Properties - Speed 


%% Comparing Wave Properties - Direction
figure();
subplot(1,2,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesHit,1500,2500));
title('Baseline');
subplot(1,2,2);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesOpto.wavesHit,1500,2500));
title('Opto');
sgtitle('Wave Direction')

[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);

% Plotting Wave Rasters colorcoded by Direction
figure();
Spikes = vertcat(selectWaves(WavesBaseline.wavesHit,1,3000).waveStart);
prop = arrayfun(@(s) s.waveDir, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
subplot(2,1,1);
rasterPlotPropColor(Spikes,prop,prop2);

Spikes = vertcat(selectWaves(WavesOpto.wavesHit,1,3000).waveStart);
prop = arrayfun(@(s) s.waveDir, selectWaves(WavesOpto.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesOpto.wavesHit,1,3000), 'UniformOutput', false);
subplot(2,1,2);
rasterPlotPropColor(Spikes,prop,prop2);

%% Plotting Wave Rasters

%  Waves during hits 
plotWaveRaster(WavesBaseline.wavesHit,WavesOpto.wavesHit,IntanBehaviourBaseline.cueHitTrace,IntanBehaviourOpto.cueHitTrace,parameters);
sgtitle('Wave Rasters for Hits: Baseline vs Opto');
subplot(4,1,1);
RTTraceTime = (IntanBehaviourBaseline.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourBaseline.cueHitTrace,2),'.r');
subplot(4,1,3);
RTTraceTime = (IntanBehaviourOpto.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourOpto.cueHitTrace,2),'.r');

% Waves during misses
plotWaveRaster(WavesBaseline.wavesMiss,WavesOpto.wavesMiss,IntanBehaviourBaseline.cueMissTrace,IntanBehaviourOpto.cueMissTrace,parameters);

% Waves during FA
plotWaveRaster(WavesBaseline.wavesFA,WavesOpto.wavesFA,IntanBehaviourBaseline.missTrace,IntanBehaviourOpto.missTrace,parameters);

% Waves during MIFA
plotWaveRaster(WavesBaseline.wavesMIFA,WavesOpto.wavesMIFA,IntanBehaviourBaseline.MIFATrace,IntanBehaviourOpto.MIFATrace,parameters);


