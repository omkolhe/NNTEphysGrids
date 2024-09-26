%% Plotting temp 
figure,
plot(IntanBehaviour.time/60,IntanBehaviour.tempTrace);


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
ylabel('Trials');xlabel('Time (in ms)');xlim([0 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);%ylim([1 20]);
set(gca,'TickDir','out','fontsize',14'); box off;

subplot(4,1,2)
bar((sum(wavesHitPresent,1)/size(IntanBehaviour.cueHitTrace,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
xline((mean(IntanBehaviour.reactionTime,'all')*parameters.Fs + parameters.windowBeforeCue*parameters.Fs+1),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylim([0 0.4]);xlim([0 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);
ylabel('Wave probability');xlabel('Time (in ms)');set(gca,'TickDir','out','fontsize',14'); box off;

subplot(4,1,3);
title('Waves During Miss Trials')
rasterPlotColor(wavesMissPresent,tempMiss);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([0 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);%ylim([1 20]);
set(gca,'TickDir','out','fontsize',14'); box off;

subplot(4,1,4)
bar((sum(wavesMissPresent,1)/size(IntanBehaviour.cueMissTrace,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
ylim([0 0.4]);xlim([0 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs])
ylabel('Wave probability');xlabel('Time (in ms)')
set(gca,'TickDir','out','fontsize',14'); box off;

%% Hit waves comparison 
baselineSpeed = [];
cooledSpeed = [];
recoverySpeed = [];
recoveryflag = 0;
waveStruct = Waves.wavesHit;
for i = 1:size(waveStruct,2)
    if IntanBehaviour.cueHitTrace(i).temp <= 12
        cooledSpeed = [cooledSpeed,waveStruct(i).speed];
        recoveryflag =1;
    elseif IntanBehaviour.cueHitTrace(i).temp >= 20 && recoveryflag == 0
        baselineSpeed = [baselineSpeed,waveStruct(i).speed];
    elseif IntanBehaviour.cueHitTrace(i).temp >= 15 && recoveryflag == 1
        recoverySpeed = [recoverySpeed,waveStruct(i).speed];
    end
end

ranksum(baselineSpeed',cooledSpeed')
ranksum(cooledSpeed',recoverySpeed')
ranksum(baselineSpeed',recoverySpeed')

temp = [baselineSpeed';cooledSpeed';recoverySpeed'];
templabels = cellstr([repmat('Baseline',size(baselineSpeed,2),1);repmat('Cooled  ',size(cooledSpeed,2),1);repmat('Recovery',size(recoverySpeed,2),1)]);
colors = [166/255 14/255 90/255;53/255 189/255 206/255;0.8500 0.3250 0.0980];
figure,violinplot(temp,templabels,'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
ylabel('Wave Speed (cm/s)'); title('Cooling M2');
% xtix = {'Baseline','Cooled','Recovery'}; xtixloc = [1 2 3]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);set(gca,'TickDir','out','fontsize',14');
set(gca,'TickDir','out','fontsize',14');box off;
%% Plotting Speed histogram

h1 = histfit(baselineSpeed,80,'lognormal');
h2 = histfit(cooledSpeed,80,'lognormal');
figure('Name','Histogram of wave speeds: Baseline vs Cooled');
bar(h1(1).XData, h1(1).YData/max(h1(1).YData),'BarWidth',2,'FaceColor',[0.7490 0 0.3490],'FaceAlpha',0.6,'EdgeColor','none');hold on;
plot(h1(2).XData, h1(2).YData/max(h1(1).YData), 'r', 'LineWidth',2)
xline(mean(baselineSpeed,'all'),'-r',{num2str(mean(baselineSpeed,'all')) ' cm/s'});
xlabel('Wave speed in cm/s');title('Wave Speed');box off;ylabel('Normalized Frequency');
bar(h2(1).XData, h2(1).YData/max(h2(1).YData),'BarWidth',2,'FaceColor',[0 1 0.9804],'FaceAlpha',0.8,'EdgeColor','none');hold on;
plot(h2(2).XData, h2(2).YData/max(h2(1).YData), 'b', 'LineWidth',2)
xline(mean(cooledSpeed,'all'),'-r',{num2str(mean(cooledSpeed,'all')) ' cm/s'});ylim([0 1.2])
set(gca,'TickDir','out','fontsize',14');


%%
RTbaseline = cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour1.cueHitTrace, 'UniformOutput', false));
RTcooling = cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace, 'UniformOutput', false));

ranksum(RTcooling,RTbaseline)
data = RTbaseline';
RTcooling';

figure,customBoxplot(data);
ylabel('Reaction Time (s)');


baselineSpeed = [];
cooledSpeed = [];

for i = 1:size(Waves1.wavesHit,2)
    baselineSpeed = [baselineSpeed,Waves1.wavesHit(i).speed];
end

for i = 1:size(Waves.wavesHit,2)
    cooledSpeed = [cooledSpeed,Waves.wavesHit(i).speed];
end

ranksum(baselineSpeed,cooledSpeed)
ttest2(baselineSpeed,cooledSpeed)
data = baselineSpeed';
cooledSpeed';
figure,customBoxplot(data);

median(baselineSpeed)-median(cooledSpeed)
