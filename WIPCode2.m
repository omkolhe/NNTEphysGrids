%% PPL across grids 

phaseResetPPL = zeros(size(IntanBehaviour.cueHitTrace,2),size(IntanBehaviour.cueHitTrace(1).xgp,3));
phaseResetPA = zeros(size(IntanBehaviour.cueHitTrace,2),size(IntanBehaviour.cueHitTrace(1).xgp,3));

% N = ceil(exp(0.626 + (0.4*log((parameters.rows*parameters.cols)-1))));
N = 36;
Hmax = log2(N);

for i=1:size(IntanBehaviour.cueHitTrace,2)
    phi = reshape(angle(IntanBehaviour.cueHitTrace(i).xgp),parameters.cols*parameters.rows,[]);
    phi = removeNaNRows(phi);
    for j=1:size(IntanBehaviour.cueHitTrace(1).xgp,3)   
        phaseResetPPL(i,j) = 100*(1-(getShannonEntropy(squeeze(phi(:,j))',N,-pi,pi)/Hmax));
        phaseResetPA(i,j) = size(phi,1) - abs(sum(exp(1i*phi(:,j))));
    end
end


figure();
subplot(2,1,1);
imagesc(IntanBehaviour.cueHitTrace(1).time,1:size(IntanBehaviour.cueHitTrace,2),phaseResetPPL);
colormap(hot);axis xy;
subplot(2,1,2);
imagesc(IntanBehaviour.cueHitTrace(1).time,1:size(IntanBehaviour.cueHitTrace,2),phaseResetPA);
colormap(hot);axis xy;

%%
trialno = 6;
figure();
% plot(IntanBehaviour.cueHitTrace(trialno).time,phaseResetPPL(trialno,:));
hold on;xline(mean(IntanBehaviour.reactionTime(trialno)));
plot(IntanBehaviour.cueHitTrace(trialno).time,phaseResetPA(trialno,:));
% yyaxis right; plot(IntanBehaviour.cueHitTrace(trialno).time,Waves.wavesHit(1).wavePresent);

[pks,loc] = findpeaks(phaseResetPA(trialno,:));
scatter((loc-1501)/1000,pks)

postCueLocs = find(loc>parameters.windowBeforeCue*parameters.Fs);
pks = pks(postCueLocs);loc = loc(postCueLocs);
highPeaks = find(pks>=16);
pks = pks(highPeaks);loc = loc(highPeaks);


%% 
firstPeakLoc = zeros(1,size(IntanBehaviour.cueHitTrace,2));
for i=1:size(IntanBehaviour.cueHitTrace,2)
    [pks,loc] = findpeaks(phaseResetPA(i,:));
    postCueLocs = find(loc>parameters.windowBeforeCue*parameters.Fs);
    pks = pks(postCueLocs);loc = loc(postCueLocs);
    highPeaks = find(pks>=10);
    pks = pks(highPeaks);loc = loc(highPeaks);
    firstPeakLoc(i) = loc(1)-(parameters.windowBeforeCue*parameters.Fs);
end

mdl = fitlm(firstPeakLoc,IntanBehaviour.reactionTime)
figure,plot(mdl);

%%
combIntanBehaviour(1).cueHitTrace = rmfield(combIntanBehaviour(1).cueHitTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).cueMissTrace = rmfield(combIntanBehaviour(1).cueMissTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).hitTrace = rmfield(combIntanBehaviour(1).hitTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).missTrace = rmfield(combIntanBehaviour(1).missTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).MIHitTrace = rmfield(combIntanBehaviour(1).MIHitTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).MIFATrace = rmfield(combIntanBehaviour(1).MIFATrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});



combIntanBehaviour(1).cueHitTrace = rmfield(combIntanBehaviour(1).cueHitTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).cueMissTrace = rmfield(combIntanBehaviour(1).cueMissTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).hitTrace = rmfield(combIntanBehaviour(1).hitTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).missTrace = rmfield(combIntanBehaviour(1).missTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).MIHitTrace = rmfield(combIntanBehaviour(1).MIHitTrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});
combIntanBehaviour(1).MIFATrace = rmfield(combIntanBehaviour(1).MIFATrace, {'rawLFPProbe', 'xfProbe', 'xgpProbe', 'wtProbe'});


%%
xgp = arrayfun(@(s) squeeze(angle(s.xgp(5,2,:)))', IntanBehaviour.cueHitTrace, 'UniformOutput', false)';

allPhase = zeros(size(xgp,1),size(xgp{1,1},2));
for i=1:size(xgp,1)
    allPhase(i,:) = xgp{i,1};
end

figure();
imagesc(allPhase);
axis xy;

figure();
plot(allPhase(1,:));

trialno = 15;
figure();
subplot(4,1,1);
stack_plot(reshape(angle(IntanBehaviour.cueHitTrace(trialno).xgp),[parameters.rows*parameters.cols],[]),0,10,parameters.Fs);
hold on; xline(1.5);xline(IntanBehaviour.cueHitTrace(trialno).reactionTime+1.5,'-r');
subplot(4,1,2);
stack_plot(reshape(angle(IntanBehaviour.cueHitTrace(trialno+1).xgp),[parameters.rows*parameters.cols],[]),0,10,parameters.Fs);
hold on; xline(1.5);xline(IntanBehaviour.cueHitTrace(trialno+1).reactionTime+1.5,'-r');
subplot(4,1,3);
stack_plot(reshape(angle(IntanBehaviour.cueHitTrace(trialno+2).xgp),[parameters.rows*parameters.cols],[]),0,10,parameters.Fs);
hold on; xline(1.5);xline(IntanBehaviour.cueHitTrace(trialno+2).reactionTime+1.5,'-r');
subplot(4,1,4);
stack_plot(reshape(angle(IntanBehaviour.cueHitTrace(trialno+3).xgp),[parameters.rows*parameters.cols],[]),0,10,parameters.Fs);
hold on; xline(1.5);xline(IntanBehaviour.cueHitTrace(trialno+3).reactionTime+1.5,'-r');


%% 

PAHitz = zscore(PPLGrid.PPLMIHit,0,[1,2,3]);
PAMissz = zscore(PPLGrid.PPLMIFA,0,[1,2,3]);

figure();
title("Phase Alignment averaged across Electrodes")
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAHitz,[1 2])),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueMissTrace(1).time,squeeze(nanmean(PAMissz,[1 2])),'-k','LineWidth',1);
% plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAFA,[1 2])),'-k','LineWidth',1);
ylabel("Phase Alignment"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Phase Alignment for Hits');box off;legend('Hits','Miss');

PPLHitz = zscore(PPLGrid.PPLHit,0,[1,2,3]);
PPLMissz = zscore(PPLGrid.PPLMiss,0,[1,2,3]);

figure();
title("Phase Alignment averaged across Electrodes")
plot(IntanBehaviour.cueHitTrace(1).time,smooth(squeeze(nanmean(PPLHitz,[1 2])),50,'sgolay',5),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueMissTrace(1).time,smooth(squeeze(nanmean(PPLMissz,[1 2])),50,'sgolay',5),'-k','LineWidth',1);
% plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAFA,[1 2])),'-k','LineWidth',1);
ylabel("z-score"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m',' RT','LabelVerticalAlignment','top');
title('Phase Locking for Hits vs Misses');box off;legend('Hits','Misses');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14'); ylim([-1.2 5]);

PPLHitz = zscore(PPLGrid.PPLMIHit,0,[1,2,3]);
PPLMissz = zscore(PPLGrid.PPLMIFA,0,[1,2,3]);

figure();
title("Phase Alignment averaged across Electrodes")
plot(IntanBehaviour.cueHitTrace(1).time,smooth(squeeze(nanmean(PPLHitz,[1 2])),50,'sgolay',5),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueMissTrace(1).time,smooth(squeeze(nanmean(PPLMissz,[1 2])),50,'sgolay',5),'-k','LineWidth',1);
% plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAFA,[1 2])),'-k','LineWidth',1);
ylabel("z-score"); xlabel("Time (s)");
xline(0,'--r','MI','LabelVerticalAlignment','top');
xline(-mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue','LabelVerticalAlignment','top');
title('Phase Locking for Hits vs FAs');box off;legend('Hits','FAs');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');ylim([-1.2 5]);



%%

wavesHitPresent = vertcat(Waves.wavesHit.wavePresent);
wavesMissPresent = vertcat(Waves.wavesMiss.wavePresent);
wavesHitStart = vertcat(Waves.wavesHit.waveStart);
wavesMissStart = vertcat(Waves.wavesMiss.waveStart);

RTTraceTime = IntanBehaviour.reactionTime;
kernal = [0 1 1 0;0 1 1 0];

figure();
subplot(4,1,1);
title('Waves During Hit Trials')
% rasterPlot(wavesHitPresent);
imagesc(IntanBehaviour.cueHitTrace(1).time,1:size(wavesHitPresent),conv2(wavesHitPresent,kernal,'valid'));colormap(flip(gray));axis xy;hold on;
xline(0,'--r','Cue','LabelVerticalAlignment','top');caxis([0 1]);
% plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','RT','LabelVerticalAlignment','top');
% xline(mean(IntanBehaviour.rewardTime,'all'),'--m','Reward','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([-0.5 1.5]);box off; ylim([1 100]);set(gca,'TickDir','out','fontsize',14')
subplot(4,1,2)
bar(IntanBehaviour.cueHitTrace(1).time,(sum(wavesHitPresent,1)/size(IntanBehaviour.cueHitTrace,2)));
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','RT','LabelVerticalAlignment','top');ylim([0 0.3]);
% xline(mean(IntanBehaviour.rewardTime,'all'),'--m','Reward','LabelVerticalAlignment','top');
ylabel('Wave probability');xlabel('Time (in ms)');xlim([-0.5 1.5]);box off;set(gca,'TickDir','out','fontsize',14')
subplot(4,1,3);
title('Waves During Miss Trials')
% rasterPlot(wavesMissPresent);
imagesc(IntanBehaviour.cueMissTrace(1).time,1:size(wavesMissPresent),conv2(wavesMissPresent,kernal,'valid'));colormap(flip(gray));axis xy;hold on;
xline(0,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([-0.5 1.5]);box off;set(gca,'TickDir','out','fontsize',14');ylim([1 80]);caxis([0 1]);
subplot(4,1,4)
bar(IntanBehaviour.cueMissTrace(1).time,(sum(wavesMissPresent,1)/size(IntanBehaviour.cueMissTrace,2)));
xline(0,'--r','Cue','LabelVerticalAlignment','top');
ylim([0 0.3]);xlim([-0.5 1.5]);box off;
ylabel('Wave probability');xlabel('Time (in ms)')
set(gca,'TickDir','out','fontsize',14')

%%

wavesHitPresent = vertcat(Waves.wavesMIHit.wavePresent);
wavesMissPresent = vertcat(Waves.wavesMIFA.wavePresent);
wavesHitStart = vertcat(Waves.wavesMIHit.waveStart);
wavesMissStart = vertcat(Waves.wavesMIFA.waveStart);

RTTraceTime = IntanBehaviour.reactionTime;
kernal = [0 1 1 0;0 1 1 0];

figure();
subplot(4,1,1);
title('Waves During Hit Trials')
% rasterPlot(wavesHitPresent);
imagesc(IntanBehaviour.MIHitTrace(1).time,1:size(wavesHitPresent),conv2(wavesHitPresent,kernal,'valid'));colormap(flip(gray));axis xy;hold on;
xline(0,'--r','MI','LabelVerticalAlignment','top');caxis([0 1]);
% plot(RTTraceTime,1:size(IntanBehaviour.cueHitTrace,2),'.r');
% xline(mean(IntanBehaviour.reactionTime,'all'),'--m','RT','LabelVerticalAlignment','top');
% xline(mean(IntanBehaviour.rewardTime,'all'),'--m','Reward','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([-0.5 1.5]);box off; ylim([1 19]);set(gca,'TickDir','out','fontsize',14')
subplot(4,1,2)
bar(IntanBehaviour.MIHitTrace(1).time,(sum(wavesHitPresent,1)/size(IntanBehaviour.MIHitTrace,2)));
xline(0,'--r','MI','LabelVerticalAlignment','top');
% xline(mean(IntanBehaviour.reactionTime,'all'),'--m','RT','LabelVerticalAlignment','top');ylim([0 0.3]);
% xline(mean(IntanBehaviour.rewardTime,'all'),'--m','Reward','LabelVerticalAlignment','top');
ylabel('Wave probability');xlabel('Time (in ms)');xlim([-0.5 1.5]);box off;set(gca,'TickDir','out','fontsize',14')
subplot(4,1,3);
title('Waves During Miss Trials')
% rasterPlot(wavesMissPresent);
imagesc(IntanBehaviour.MIFATrace(1).time,1:size(wavesMissPresent),conv2(wavesMissPresent,kernal,'valid'));colormap(flip(gray));axis xy;hold on;
xline(0,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([-0.5 1.5]);box off;set(gca,'TickDir','out','fontsize',14');ylim([1 80]);caxis([0 1]);
subplot(4,1,4)
bar(IntanBehaviour.MIFATrace(1).time,(sum(wavesMissPresent,1)/size(IntanBehaviour.MIFATrace,2)));
xline(0,'--r','Cue','LabelVerticalAlignment','top');
ylim([0 0.3]);xlim([-0.5 1.5]);box off;
ylabel('Wave probability');xlabel('Time (in ms)')
set(gca,'TickDir','out','fontsize',14')


%%
wavesHitPresent = vertcat(Waves1.wavesHit.wavePresent);
wavesMissPresent = vertcat(Waves1.wavesMiss.wavePresent);
wavesHitStart = vertcat(Waves1.wavesHit.waveStart);
wavesMissStart = vertcat(Waves1.wavesMiss.waveStart);

RTTraceTime = (IntanBehaviour1.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);

figure();
subplot(2,1,1);
title('Waves During Hit Trials')
rasterPlot(wavesHitPresent);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
plot(RTTraceTime,1:size(IntanBehaviour1.cueHitTrace,2),'.r');
xline((mean(IntanBehaviour1.reactionTime,'all')*parameters.Fs + parameters.windowBeforeCue*parameters.Fs+1),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);ylim([1 20]);
set(gca,'TickDir','out','fontsize',14'); box off;
subplot(2,1,2)
bar((sum(wavesHitPresent(1:19,:),1)/19));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
xline((mean(IntanBehaviour1.reactionTime,'all')*parameters.Fs + parameters.windowBeforeCue*parameters.Fs+1),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylim([0 0.3]);xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);
ylabel('Wave probability');xlabel('Time (in ms)');set(gca,'TickDir','out','fontsize',14'); box off;


%%
avgSpeedBaseline = mean(WavesBaseline.speed,'all');
figure('Name','Histogram of wave speeds');
subplot(2,1,1)
h1 = histfit(WavesBaseline.speed,100,'kernel');
xline(avgSpeedBaseline,'-r',{'Mean speed = ' num2str(avgSpeedBaseline) ' cm/s'});
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed - Baseline');

avgSpeedCool = mean(WavesCool.speed,'all');
subplot(2,1,2)
h2 = histfit(WavesCool.speed,100,'kernel');
xline(avgSpeedCool,'-r',{'Mean speed = ' num2str(avgSpeedCool) ' cm/s'});
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed - Cooled');


ranksum(log(WavesCool.speed), log(WavesBaseline.speed))
