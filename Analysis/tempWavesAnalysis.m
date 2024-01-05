%% Plotting temp 
figure,
plot(IntanBehaviour.time,IntanBehaviour.tempTrace);


%% Waves raster plot - color coded as temperature 
wavesHitPresent = vertcat(Waves.wavesHit.wavePresent);
wavesMissPresent = vertcat(Waves.wavesMiss.wavePresent);
wavesHitStart = vertcat(Waves.wavesHit.waveStart);
wavesMissStart = vertcat(Waves.wavesMiss.waveStart);
tempHit = vertcat(IntanBehaviour.cueHitTrace.temp);
tempMiss = vertcat(IntanBehaviour.cueMissTrace.temp);

RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);

figure();
subplot(4,1,1);
title('Waves During Hit Trials')
rasterPlotColor(wavesHitPresent,tempHit);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');
xline((mean(IntanBehaviour.reactionTime,'all')*parameters.Fs + parameters.windowBeforeCue*parameters.Fs+1),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);%ylim([1 20]);
set(gca,'TickDir','out','fontsize',14'); box off;

subplot(4,1,2)
bar((sum(wavesHitPresent,1)/size(IntanBehaviour.cueHitTrace,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
xline((mean(IntanBehaviour.reactionTime,'all')*parameters.Fs + parameters.windowBeforeCue*parameters.Fs+1),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylim([0 0.3]);xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);
ylabel('Wave probability');xlabel('Time (in ms)');set(gca,'TickDir','out','fontsize',14'); box off;

subplot(4,1,3);
title('Waves During Miss Trials')
rasterPlotColor(wavesMissPresent,tempMiss);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);%ylim([1 20]);
set(gca,'TickDir','out','fontsize',14'); box off;

subplot(4,1,4)
bar((sum(wavesMissPresent,1)/size(IntanBehaviour.cueMissTrace,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
ylim([0 0.3]);xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs])
ylabel('Wave probability');xlabel('Time (in ms)')
set(gca,'TickDir','out','fontsize',14'); box off;

%% Hit waves comparison 
baselineSpeed = [];
cooledSpeed = [];

for i = 1:size(Waves.wavesHit,2)
    if IntanBehaviour.cueHitTrace(i).temp <= 20
        cooledSpeed = [cooledSpeed,Waves.wavesHit(i).speed];
    elseif IntanBehaviour.cueHitTrace(i).temp >= 25
        baselineSpeed = [baselineSpeed,Waves.wavesHit(i).speed];
    end
end

ranksum(baselineSpeed,cooledSpeed)
data = baselineSpeed';

figure,customBoxplot(data);
