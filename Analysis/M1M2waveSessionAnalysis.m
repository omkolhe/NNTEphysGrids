%% Plotting rasters 
%  Waves during hits for M1 and M2
plotWaveRaster(WavesM1.wavesHit,WavesM2.wavesHit,IntanBehaviourM1.cueHitTrace,IntanBehaviourM2.cueHitTrace,parameters);
sgtitle('Wave Rasters for Hits in M1 and M2');
subplot(4,1,1);
RTTraceTime = (IntanBehaviourM1.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourM1.cueHitTrace,2),'.r');
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
xline(median(RTTraceTime),'--r','RT','LabelVerticalAlignment','top');
subplot(4,1,2);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
subplot(4,1,3);
RTTraceTime = (IntanBehaviourM2.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviourM2.cueHitTrace,2),'.r');
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
subplot(4,1,4);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
xlim([1000 3000]);

%  Waves during misses for M1 and M2
plotWaveRaster(WavesM1.wavesMiss,WavesM2.wavesMiss,IntanBehaviourM1.cueMissTrace,IntanBehaviourM2.cueMissTrace,parameters);
sgtitle('Wave Rasters for Misses in M1 and M2');
subplot(4,1,1);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
subplot(4,1,2);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
subplot(4,1,3);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
subplot(4,1,4);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
xlim([1000 3000]);

% Waves in M1 and M2 during MI Hits
plotWaveRaster(WavesM1.wavesMIHit,WavesM2.wavesMIHit,IntanBehaviourM1.MIHitTrace,IntanBehaviourM2.MIHitTrace,parameters);
sgtitle('Wave Rasters for MIHits in M1 and M2');
subplot(4,1,1);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
subplot(4,1,2);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
subplot(4,1,3);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
subplot(4,1,4);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
xlim([1000 3000]);

% Waves in M1 and M2 during MI FAs
plotWaveRaster(WavesM1.wavesMIFA,WavesM2.wavesMIFA,IntanBehaviourM1.MIFATrace,IntanBehaviourM2.MIFATrace,parameters);
sgtitle('Wave Rasters for MIFAs in M1 and M2');
subplot(4,1,1);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
subplot(4,1,2);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
subplot(4,1,3);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
subplot(4,1,4);
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
xlim([1000 3000]);


%% Wave Direction
% Waves Hits vs Miss
W = vertcat(selectWaves(WavesM1.wavesHit,1,3000).waveStart);figure();
prop = arrayfun(@(s) s.waveDir, selectWaves(WavesM1.wavesHit,1,3000), 'UniformOutput', false);
ax1 = subplot(2,1,1);
rasterPlotPropColor(W,prop,[],1);
RTTraceTime = (IntanBehaviourM1.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');xlim([1 size(W,2)]);
title('Wave Hits - M1')

W = vertcat(selectWaves(WavesM2.wavesHit,1,3000).waveStart);
prop = arrayfun(@(s) s.waveDir, selectWaves(WavesM2.wavesHit,1,3000), 'UniformOutput', false);
ax2 = subplot(2,1,2);
rasterPlotPropColor(W,prop,[],1);
RTTraceTime = (IntanBehaviourM1.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');xlim([1 size(W,2)]);
title('Wave Hits - M2')
linkaxes([ax1,ax2],'x');

% Plotting Wave Rasters colorcoded by Direction
srt = 1500;
stp = 2200;
figure();
subplot(1,2,1);
[dirCombHit,~] = plotWaveDirection(selectWaves(WavesM1.wavesHit(1:end),srt,stp),36,[]);
title('Hits');
subplot(1,2,2);
[dirCombMiss,~] = plotWaveDirection(selectWaves(WavesM1.wavesMiss(1:end),srt,stp),36,[]);
title('Miss');
sgtitle('Wave Direction')

[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);

% Hits
srt = 1;stp = 1500;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesHit,srt,stp),36,[]);
title('M1: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesHit,srt,stp),36,[]);
title('M2: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 2200;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesHit,srt,stp),36,[]);
title('M1: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesHit,srt,stp),36,[]);
title('M2: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 2201;stp = 3001;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesHit,srt,stp),36,[]);
title('M1: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesHit,srt,stp),36,[]);
title('M2: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction')

% Miss
srt = 1;stp = 1500;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesMiss,srt,stp),36,[]);
title('M1: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesMiss,srt,stp),36,[]);
title('M2: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 2200;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesMiss,srt,stp),36,[]);
title('M1: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesMiss,srt,stp),36,[]);
title('M2: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 2201;stp = 3001;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesMiss,srt,stp),36,[]);
title('M1: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesMiss,srt,stp),36,[]);
title('M2: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction - Miss')

 % MI Hits
srt = 1;stp = 1500;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesMIHit,srt,stp),36,[]);
title('M1: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesMIHit,srt,stp),36,[]);
title('M2: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 2200;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesMIHit,srt,stp),36,[]);
title('M1: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesMIHit,srt,stp),36,[]);
title('M2: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 2201;stp = 3001;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesMIHit,srt,stp),36,[]);
title('M1: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesMIHit,srt,stp),36,[]);
title('M2: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction - MIHit')

% MI FAs
srt = 1;stp = 1500;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesMIFA,srt,stp),36,[]);
title('M1: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesMIFA,srt,stp),36,[]);
title('M2: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 2200;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesMIFA,srt,stp),36,[]);
title('M1: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesMIFA,srt,stp),36,[]);
title('M2: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 2201;stp = 3001;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(WavesM1.wavesMIFA,srt,stp),36,[]);
title('M1: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombOpto,~] = plotWaveDirection(selectWaves(WavesM2.wavesMIFA,srt,stp),36,[]);
title('M2: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction - MIFA')

%% Wave speed 
% Plotting wave speed as function of time - M1
nPoints = 30; interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
waveAvgFreq = zeros(4,nPoints);
for i=1:nPoints
    st = (i-1)*interval + 1;
    sp = (i)*interval + 1;
    WaveSpeed(i).speedHit = horzcat(selectWaves(WavesM1.wavesHit,st,sp).speed);
    WaveSpeed(i).speedMiss = horzcat(selectWaves(WavesM1.wavesMiss,st,sp).speed);
    WaveSpeed(i).speedMIHit = horzcat(selectWaves(WavesM1.wavesMIHit,st,sp).speed);
    WaveSpeed(i).speedMIFA = horzcat(selectWaves(WavesM1.wavesMIFA,st,sp).speed);
end

t = interval:interval:interval*nPoints;
y = cell2mat(arrayfun(@(s) mean(s.speedHit,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedHit,0,'all','omitnan')/sqrt(numel(s.speedHit)),WaveSpeed,'UniformOutput',false));
figure;hold on;
h1 = errorbar(t,y,err,'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
y = cell2mat(arrayfun(@(s) mean(s.speedMiss,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMiss,0,'all','omitnan')/sqrt(numel(s.speedMiss)),WaveSpeed,'UniformOutput',false));
h2 = errorbar(t,y,err,'Color', [0.7 0.7 0.7],'LineWidth',1.5);
y = cell2mat(arrayfun(@(s) mean(s.speedMIHit,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMIHit,0,'all','omitnan')/sqrt(numel(s.speedMIHit)),WaveSpeed,'UniformOutput',false));
h3 = errorbar(t,y,err,'Color', [0 0 1 0.4],'LineWidth',1.5);
y = cell2mat(arrayfun(@(s) mean(s.speedMIFA,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMIFA,0,'all','omitnan')/sqrt(numel(s.speedMIFA)),WaveSpeed,'UniformOutput',false));
h4 = errorbar(t,y,err,'Color', [1 0 0 0.4],'LineWidth',1.5);

xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Average Wave Frequency (Hz)');
legend([h1 h2 h3 h4],'Hits','Misses','MIHits','MIFAs','Location','best'); %ylim([5 15]);
title('Wave Speed - M1')

% Plotting wave speed as function of time - M2
nPoints = 30; interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
waveAvgFreq = zeros(4,nPoints);
for i=1:nPoints
    st = (i-1)*interval + 1;
    sp = (i)*interval + 1;
    WaveSpeed(i).speedHit = horzcat(selectWaves(WavesM2.wavesHit,st,sp).speed);
    WaveSpeed(i).speedMiss = horzcat(selectWaves(WavesM2.wavesMiss,st,sp).speed);
    WaveSpeed(i).speedMIHit = horzcat(selectWaves(WavesM2.wavesMIHit,st,sp).speed);
    WaveSpeed(i).speedMIFA = horzcat(selectWaves(WavesM2.wavesMIFA,st,sp).speed);
end

t = interval:interval:interval*nPoints;
y = cell2mat(arrayfun(@(s) mean(s.speedHit,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedHit,0,'all','omitnan')/sqrt(numel(s.speedHit)),WaveSpeed,'UniformOutput',false));
figure;hold on;
h1 = errorbar(t,y,err,'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
y = cell2mat(arrayfun(@(s) mean(s.speedMiss,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMiss,0,'all','omitnan')/sqrt(numel(s.speedMiss)),WaveSpeed,'UniformOutput',false));
h2 = errorbar(t,y,err,'Color', [0.7 0.7 0.7],'LineWidth',1.5);
y = cell2mat(arrayfun(@(s) mean(s.speedMIHit,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMIHit,0,'all','omitnan')/sqrt(numel(s.speedMIHit)),WaveSpeed,'UniformOutput',false));
h3 = errorbar(t,y,err,'Color', [0 0 1 0.4],'LineWidth',1.5);
y = cell2mat(arrayfun(@(s) mean(s.speedMIFA,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMIFA,0,'all','omitnan')/sqrt(numel(s.speedMIFA)),WaveSpeed,'UniformOutput',false));
h4 = errorbar(t,y,err,'Color', [1 0 0 0.4],'LineWidth',1.5);

xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Average Wave Frequency (Hz)');
legend([h1 h2 h3 h4],'Hits','Misses','MIHits','MIFAs','Location','best'); %ylim([5 15]);
title('Wave Speed - M2')

%% Wave PDG 

% Waves PGD in M1
PGDHits = vertcat(WavesM1.wavesHit.PGD);
PGDMiss = vertcat(WavesM1.wavesMiss.PGD);
PGDMIHit = vertcat(WavesM1.wavesMIHit.PGD);
PGDMIFA = vertcat(WavesM1.wavesMIFA.PGD);

figure();
subplot(2,1,1);
hold on;
plot(IntanBehaviourM1.cueHitTrace(1).time,smooth(mean(PGDHits,1),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
plot(IntanBehaviourM1.cueHitTrace(1).time,smooth(mean(PGDHits,1)-(std(PGDHits,0,1)/sqrt(size(PGDHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviourM1.cueHitTrace(1).time,smooth(mean(PGDHits,1)+(std(PGDHits,0,1)/sqrt(size(PGDHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);

plot(IntanBehaviourM1.cueMissTrace(1).time,smooth(mean(PGDMiss,1),100,'sgolay',4),'Color', [0.7 0.7 0.7],'LineWidth',1.5);
plot(IntanBehaviourM1.cueMissTrace(1).time,smooth(mean(PGDMiss,1)-(std(PGDMiss,0,1)/sqrt(size(PGDMiss,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
plot(IntanBehaviourM1.cueMissTrace(1).time,smooth(mean(PGDMiss,1)+(std(PGDMiss,0,1)/sqrt(size(PGDMiss,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);

ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviourM1.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Misses');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');

subplot(2,1,2);
hold on;
plot(IntanBehaviourM1.MIHitTrace(1).time,smooth(mean(PGDMIHit,1),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
plot(IntanBehaviourM1.MIHitTrace(1).time,smooth(mean(PGDMIHit,1)-(std(PGDMIHit,0,1)/sqrt(size(PGDMIHit,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviourM1.MIHitTrace(1).time,smooth(mean(PGDMIHit,1)+(std(PGDMIHit,0,1)/sqrt(size(PGDMIHit,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);

plot(IntanBehaviourM1.MIFATrace(1).time,smooth(mean(PGDMIFA,1),100,'sgolay',4),'Color', [0.7 0.7 0.7],'LineWidth',1.5);
plot(IntanBehaviourM1.MIFATrace(1).time,smooth(mean(PGDMIFA,1)-(std(PGDMIFA,0,1)/sqrt(size(PGDMIFA,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
plot(IntanBehaviourM1.MIFATrace(1).time,smooth(mean(PGDMIFA,1)+(std(PGDMIFA,0,1)/sqrt(size(PGDMIFA,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);

ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','MI','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('MIHits','MIFAs');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');

% Waves PGD in M2
PGDHits = vertcat(WavesM2.wavesHit.PGD);
PGDMiss = vertcat(WavesM2.wavesMiss.PGD);
PGDMIHit = vertcat(WavesM2.wavesMIHit.PGD);
PGDMIFA = vertcat(WavesM2.wavesMIFA.PGD);

figure();
subplot(2,1,1);
hold on;
plot(IntanBehaviourM2.cueHitTrace(1).time,smooth(mean(PGDHits,1),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
plot(IntanBehaviourM2.cueHitTrace(1).time,smooth(mean(PGDHits,1)-(std(PGDHits,0,1)/sqrt(size(PGDHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviourM2.cueHitTrace(1).time,smooth(mean(PGDHits,1)+(std(PGDHits,0,1)/sqrt(size(PGDHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);

plot(IntanBehaviourM2.cueMissTrace(1).time,smooth(mean(PGDMiss,1),100,'sgolay',4),'Color', [0.7 0.7 0.7],'LineWidth',1.5);
plot(IntanBehaviourM2.cueMissTrace(1).time,smooth(mean(PGDMiss,1)-(std(PGDMiss,0,1)/sqrt(size(PGDMiss,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
plot(IntanBehaviourM2.cueMissTrace(1).time,smooth(mean(PGDMiss,1)+(std(PGDMiss,0,1)/sqrt(size(PGDMiss,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);

ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviourM2.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Misses');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');

subplot(2,1,2);
hold on;
plot(IntanBehaviourM2.MIHitTrace(1).time,smooth(mean(PGDMIHit,1),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
plot(IntanBehaviourM2.MIHitTrace(1).time,smooth(mean(PGDMIHit,1)-(std(PGDMIHit,0,1)/sqrt(size(PGDMIHit,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviourM2.MIHitTrace(1).time,smooth(mean(PGDMIHit,1)+(std(PGDMIHit,0,1)/sqrt(size(PGDMIHit,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);

plot(IntanBehaviourM2.MIFATrace(1).time,smooth(mean(PGDMIFA,1),100,'sgolay',4),'Color', [0.7 0.7 0.7],'LineWidth',1.5);
plot(IntanBehaviourM2.MIFATrace(1).time,smooth(mean(PGDMIFA,1)-(std(PGDMIFA,0,1)/sqrt(size(PGDMIFA,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
plot(IntanBehaviourM2.MIFATrace(1).time,smooth(mean(PGDMIFA,1)+(std(PGDMIFA,0,1)/sqrt(size(PGDMIFA,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);

ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','MI','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('MIHits','MIFAs');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');
