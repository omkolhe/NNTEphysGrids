%% Plotting rasters 
%  Waves during hits vs miss
plotWaveRaster(Waves.wavesHit(1:end),Waves.wavesMiss,IntanBehaviour.cueHitTrace,IntanBehaviour.cueMissTrace,parameters);
sgtitle('Wave Rasters for Hits vs Misses');
subplot(4,1,1);
RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');

%  Waves during Hits vs FA
plotWaveRaster(Waves.wavesHitReward,Waves.wavesFA,IntanBehaviour.hitTrace,IntanBehaviour.missTrace,parameters);
sgtitle('Wave Rasters for Hits vs FA');

%  Waves during MIHits vs MIFA 
plotWaveRaster(Waves.wavesMIHit,Waves.wavesMIFA,IntanBehaviour.MIHitTrace,IntanBehaviour.MIFATrace,parameters);
sgtitle('Wave Rasters for Hits vs FA');


%% Wave Direction
% Waves Hits vs Miss
W = vertcat(selectWaves(Waves.wavesHit,1,3000).waveStart);figure();
prop = arrayfun(@(s) s.waveDir, selectWaves(Waves.wavesHit,1,3000), 'UniformOutput', false);
ax1 = subplot(2,1,1);
rasterPlotPropColor(W,prop,[],1);
RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);
% plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');xlim([1 size(W,2)]);
title('Wave Hits')

W = vertcat(selectWaves(Waves.wavesMiss,1,3000).waveStart);
prop = arrayfun(@(s) s.waveDir, selectWaves(Waves.wavesMiss,1,3000), 'UniformOutput', false);
ax2 = subplot(2,1,2);
rasterPlotPropColor(W,prop,[],1);
title('Wave Miss')
linkaxes([ax1,ax2],'x');

% Plotting Wave Rasters colorcoded by Direction
srt = 1500;
stp = 1700;
figure();
subplot(1,2,1);
[dirCombHit,~] = plotWaveDirection(selectWaves(Waves.wavesHit(1:end),srt,stp),36,[]);
title('Hits');
subplot(1,2,2);
[dirCombMiss,~] = plotWaveDirection(selectWaves(Waves.wavesMiss(1:end),srt,stp),36,[]);
title('Miss');
sgtitle('Wave Direction')

[p,~,~] = circ_kuipertest(dirCombHit, dirCombMiss,60,0);
% Print the results.
disp('Wave Direction')
disp('p-value:');
disp(p);


srt = 1000;stp = 1300;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesMIHit,srt,stp),36,[]);
title('Baseline: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMIFA,srt,stp),36,[]);
title('Opto: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1300;stp = 1700;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesMIHit,srt,stp),36,[]);
title('Baseline: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMIFA,srt,stp),36,[]);
title('Opto: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1700;stp = 3000;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesMIHit,srt,stp),36,[]);
title('Baseline: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMIFA,srt,stp),36,[]);
title('Opto: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction')



srt = 1000;stp = 1500;
figure();
subplot(2,3,1);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesHit,srt,stp),36,[]);
title('Baseline: Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,4);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMiss,srt,stp),36,[]);
title('Opto: Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1500;stp = 1800;
subplot(2,3,2);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesHit,srt,stp),36,[]);
title('Baseline: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,5);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMiss,srt,stp),36,[]);
title('Opto: Cue Evoked');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
srt = 1800;stp = 2200;
subplot(2,3,3);
[dirCombBaseline,~] = plotWaveDirection(selectWaves(Waves.wavesHit,srt,stp),36,[]);
title('Baseline: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
subplot(2,3,6);
[dirCombOpto,~] = plotWaveDirection(selectWaves(Waves.wavesMiss,srt,stp),36,[]);
title('Opto: Post Cue Spontaneous');set(gca,'TickDir','out','fontsize',14');
[p,~,~] = circ_kuipertest(dirCombBaseline, dirCombOpto,60,0);disp('Wave Direction');disp('p-value:');disp(p);
sgtitle('Wave Direction')

%% Waves speed vs time 
nPoints = 30; interval = (parameters.Fs*(parameters.windowAfterCue+parameters.windowBeforeCue))/nPoints;
waveAvgFreq = zeros(4,nPoints);
for i=1:nPoints
    st = (i-1)*interval + 1;
    sp = (i)*interval + 1;
    WaveSpeed(i).speedHit = horzcat(selectWaves(Waves.wavesHit,st,sp).speed);
    WaveSpeed(i).speedMiss = horzcat(selectWaves(Waves.wavesMiss,st,sp).speed);
    WaveSpeed(i).speedMIHit = horzcat(selectWaves(Waves.wavesMIHit,st,sp).speed);
    WaveSpeed(i).speedMIFA = horzcat(selectWaves(Waves.wavesMIFA,st,sp).speed);
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

figure();
y = cell2mat(arrayfun(@(s) mean(s.speedMIHit,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMIHit,0,'all','omitnan')/sqrt(numel(s.speedMIHit)),WaveSpeed,'UniformOutput',false));
h1 = plot(t,smoothdata(y,'gaussian'),'Color', [0.8500 0.3250 0.0980],'LineWidth',2); hold on;
plot(t,smoothdata(y-err,'gaussian'),'Color', [0.8500 0.3250 0.0980],'LineWidth',1); hold on;
plot(t,smoothdata(y+err,'gaussian'),'Color', [0.8500 0.3250 0.0980],'LineWidth',1); hold on;
y = cell2mat(arrayfun(@(s) mean(s.speedMiss,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMiss,0,'all','omitnan')/sqrt(numel(s.speedMiss)),WaveSpeed,'UniformOutput',false));
h2 = plot(t,smoothdata(y,'gaussian'),'Color', [0 0 0],'LineWidth',2); hold on;
plot(t,smoothdata(y-err,'gaussian'),'Color', [0.7 0.7 0.7],'LineWidth',1); hold on;
plot(t,smoothdata(y+err,'gaussian'),'Color', [0.7 0.7 0.7],'LineWidth',1); hold on;

xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Average Wave Speed (Hz)');
legend([h1 h2],'Hits','Misses','Location','best'); %ylim([5 15]);
title('Wave Speed - M1')
xlim([1000 3000]);box off;set(gca,'TickDir','out','fontsize',14');


figure();
y = cell2mat(arrayfun(@(s) mean(s.speedHit,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedHit,0,'all','omitnan')/sqrt(numel(s.speedHit)),WaveSpeed,'UniformOutput',false));
h1 = plot(t,smoothdata(y,'gaussian'),'Color', [0.8500 0.3250 0.0980],'LineWidth',2); hold on;
plot(t,smoothdata(y-err,'gaussian'),'Color', [0.8500 0.3250 0.0980],'LineWidth',1); hold on;
plot(t,smoothdata(y+err,'gaussian'),'Color', [0.8500 0.3250 0.0980],'LineWidth',1); hold on;
y = cell2mat(arrayfun(@(s) mean(s.speedMIFA,'all','omitnan'),WaveSpeed,'UniformOutput',false));
err = cell2mat(arrayfun(@(s) std(s.speedMIFA,0,'all','omitnan')/sqrt(numel(s.speedMIFA)),WaveSpeed,'UniformOutput',false));
h2 = plot(t,smoothdata(y,'gaussian'),'Color', [0 0 0],'LineWidth',2); hold on;
plot(t,smoothdata(y-err,'gaussian'),'Color', [0.7 0.7 0.7],'LineWidth',1); hold on;
plot(t,smoothdata(y+err,'gaussian'),'Color', [0.7 0.7 0.7],'LineWidth',1); hold on;

xline(1501,'--r','MI');xlabel('Time (ms)'); ylabel('Average Wave Speed (Hz)');
legend([h1 h2],'Hits','FAs','Location','best'); %ylim([5 15]);
title('Wave Speed - M1')
xlim([1000 3000]);box off;set(gca,'TickDir','out','fontsize',14');



