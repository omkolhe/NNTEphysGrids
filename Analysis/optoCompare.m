%% Comaparing Behaviour 
[p,h] = ranksum(IntanBehaviourBaseline.reactionTime,IntanBehaviourOpto.reactionTime)
plotBox2(IntanBehaviourBaseline.reactionTime,IntanBehaviourOpto.reactionTime);
% xL=xlim;
% yL=ylim;
% text(0.995*xL(2),0.995*yL(2),['p-val = ' num2str(p)],'HorizontalAlignment','right','VerticalAlignment','top')
ylabel('Reaction Time (s)'); title('M2 -> Th Opto');subtitle(['p-val = ' num2str(p)]);
xtix = {'Baseline','Opto'}; xtixloc = [1 2]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
set(gca,'TickDir','out','fontsize',14');

%% Comparing PA 
PABaseline = getPA(IntanBehaviourBaseline,0,1,0,parameters,0);
PAOpto = getPA(IntanBehaviourOpto,0,1,0,parameters,0);
% Smooth
figure();
subplot(2,2,[1,2])
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(squeeze(mean(PABaseline.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[0.7 0.7 0.7],'LineWidth',1.5); hold on;
plot(IntanBehaviourOpto.cueHitTrace(1).time,smooth(squeeze(mean(PAOpto.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
ylabel("Phase Alignment"); xlabel("Time (s)");
xline(0,'--k','Cue','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xline((PABaseline.PAPeakHit/parameters.Fs-parameters.windowBeforeCue),'--b','Peak','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xlim([-0.5 1.5]);box off;legend('Baseline','Opto');legend('boxoff');set(gca,'TickDir','out','fontsize',14');
title("Phase Alignment : M2 -> Th Opto");

subplot(2,2,3)
PABaselinePeak = PABaseline.PAPeakHit;
xgp = arrayfun(@(s) s.xgp, IntanBehaviourBaseline.cueHitTrace, 'UniformOutput', false);
PABaselineAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PABaselinePeak)),parameters.rows*parameters.cols,1), xgp,'UniformOutput',false)),[],1));
histogram(PABaselineAngles,18,'FaceAlpha',0.7,'FaceColor',[0.7 0.7 0.7],'Normalization','probability','EdgeColor','none');hold on;
xgp = arrayfun(@(s) s.xgp, IntanBehaviourOpto.cueHitTrace, 'UniformOutput', false);
PAOptoAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PABaselinePeak)),parameters.rows*parameters.cols,1), xgp,'UniformOutput',false)),[],1));
histogram(PAOptoAngles,18,'FaceAlpha',0.75,'FaceColor',[0.8500 0.3250 0.0980],'Normalization','probability','EdgeColor','none');box off;
xlabel('Peak PA Angle');ylabel('Probablitity');title('PA Angle at Peak');
subplot(2,2,4)
polarhistogram(PABaselineAngles,18,'FaceAlpha',0.7,'FaceColor',[0.7 0.7 0.7],'Normalization','probability','EdgeColor','none');hold on;
polarhistogram(PAOptoAngles,18,'FaceAlpha',0.75,'FaceColor',[0.8500 0.3250 0.0980],'Normalization','probability','EdgeColor','none');box off;
title('PA Angle at Peak');

[p,~,~] = circ_kuipertest(PABaselineAngles, PAOptoAngles,60,0);
disp(['Peak PA angle p-val = ' num2str(p)]);

% [p,h] = getContStat(squeeze(reshape(PABaseline.Hit,32,1,[])),squeeze(reshape(PAOpto.Hit,32,1,[])));


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
% Plotting Wave Rasters colorcoded by Speed
figure();
Spikes = vertcat(selectWaves(WavesBaseline.wavesHit,1,3000).waveStart);
prop = arrayfun(@(s) s.speed, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
subplot(2,1,1);
rasterPlotPropColor(Spikes,prop,prop2,0);
RTTraceTime = (IntanBehaviourBaseline.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourBaseline.cueHitTrace,2),'.r'); 
title('Wave Hits - Baseline')
Spikes = vertcat(selectWaves(WavesOpto.wavesHit,1,3000).waveStart);
prop = arrayfun(@(s) s.speed, selectWaves(WavesOpto.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesOpto.wavesHit,1,3000), 'UniformOutput', false);
subplot(2,1,2);
rasterPlotPropColor(Spikes,prop,prop2,0);
title('Waves Hits - Opto');
RTTraceTime = (IntanBehaviourOpto.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourOpto.cueHitTrace,2),'.r');
sgtitle('Wave rasters for Hits - Baseline vs Opto');

% Plotting wave speed histograms
srt = 1500;
stp = 2500;
speedCombBaseline = horzcat(selectWaves(WavesBaseline.wavesHit,srt,stp).speed);
speedCombOpto = horzcat(selectWaves(WavesOpto.wavesHit,srt,stp).speed);

figure('Name','Histogram of wave speeds in Baseline and Opto');
subplot(2,1,1);
h1 = histfit(speedCombBaseline,100,'lognormal');
h1(1).FaceAlpha=0.7; h1(1).FaceColor=[0.7 0.7 0.7]; h1(1).EdgeColor='none';
xline(mean(speedCombBaseline),'-r',{'Mean speed = ' num2str(mean(speedCombBaseline)) ' cm/s'});
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Baseline');box off;xlim([0 30]);
subplot(2,1,2);
h2 = histfit(speedCombOpto,100,'lognormal');
h2(1).FaceAlpha=0.7; h2(1).FaceColor=[0.8500 0.3250 0.0980]; h2(1).EdgeColor='none';
xline(mean(speedCombOpto),'-r',{'Mean speed = ' num2str(mean(speedCombOpto)) ' cm/s'});
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Opto');box off;xlim([0 30]);

figure('Name','Wave speeds in Baseline and Opto');
plotBox2(speedCombBaseline,speedCombOpto);box off;
set(gca,'XTickLabel',{'Baseline','Opto'});
ylabel('Wave speed in cm/s');ylim([0 30]);
% Perform the t-test.
[p, t] = ranksum(speedCombBaseline, speedCombOpto);
% Print the results.
disp('Wave Speed')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);
%% Comparing Wave Properties - Direction

% Plotting Wave Rasters colorcoded by Direction
srt = 1;
stp = 1500;
figure();
subplot(2,2,1);
[~,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesHit,srt,stp),[]);
title('Baseline: Spontaneous');
subplot(2,2,2);
[~,~] = plotWaveDirection(selectWaves(WavesOpto.wavesHit,srt,stp),[]);
title('Opto: Spontaneous');
sgtitle('Wave Direction')
srt = 1500;
stp = 3000;
subplot(2,2,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesBaseline.wavesHit,srt,stp),[]);
title('Baseline: Cue Evoked');
subplot(2,2,4);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesOpto.wavesHit,srt,stp),[]);
title('Opto: Cue Evoked');
sgtitle('Wave Direction')

[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);

Spikes = vertcat(selectWaves(WavesBaseline.wavesHit,1,3000).waveStart);figure();
prop = arrayfun(@(s) s.waveDir, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesBaseline.wavesHit,1,3000), 'UniformOutput', false);
ax1 = subplot(2,1,1);
rasterPlotPropColor(Spikes,prop,[],1);
RTTraceTime = (IntanBehaviourBaseline.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviourBaseline.cueHitTrace,2),'.r');xlim([1 size(Spikes,2)]);
title('Wave Hits - Baseline')

Spikes = vertcat(selectWaves(WavesOpto.wavesHit,1,3000).waveStart);
prop = arrayfun(@(s) s.waveDir, selectWaves(WavesOpto.wavesHit,1,3000), 'UniformOutput', false);
prop2 = arrayfun(@(s) s.waveDuration, selectWaves(WavesOpto.wavesHit,1,3000), 'UniformOutput', false);
ax2 = subplot(2,1,2);
rasterPlotPropColor(Spikes,prop,[],1);
title('Waves Hits - Opto');
RTTraceTime = (IntanBehaviourOpto.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviourOpto.cueHitTrace,2),'.r'); xlim([1 size(Spikes,2)]);
sgtitle('Wave rasters for Hits - Baseline vs Opto');
linkaxes([ax1,ax2],'x');
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


