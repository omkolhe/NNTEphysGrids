%% Work In Progress code - Scratch work for new fucntions

%% Trying to see if power in different bands modulate wave properties 

%% Beta power 
Beta.hitPower = arrayfun(@(s) (abs(s.xgpbeta)).^2,IntanBehaviour.cueHitTrace,'UniformOutput',false);
Beta.hitPhase = arrayfun(@(s) angle(s.xgpbeta),IntanBehaviour.cueHitTrace,'UniformOutput',false);
Gamma.hitPower = arrayfun(@(s) (abs(s.xgpgamma)).^2,IntanBehaviour.cueHitTrace,'UniformOutput',false);
Gamma.hitPhase = arrayfun(@(s) angle(s.xgpgamma),IntanBehaviour.cueHitTrace,'UniformOutput',false);
Beta.missPower = arrayfun(@(s) (abs(s.xgpbeta)).^2,IntanBehaviour.cueMissTrace,'UniformOutput',false);
Beta.missPhase = arrayfun(@(s) angle(s.xgpbeta),IntanBehaviour.cueMissTrace,'UniformOutput',false);
Gamma.missPower = arrayfun(@(s) (abs(s.xgpgamma)).^2,IntanBehaviour.cueMissTrace,'UniformOutput',false);
Gamma.missPhase = arrayfun(@(s) angle(s.xgpgamma),IntanBehaviour.cueMissTrace,'UniformOutput',false);
Beta.hitAvgPower = mean(cat(4,Beta.hitPower{:}),4);
Beta.hitAvgPhase = mean(cat(4,Beta.hitPhase{:}),4);
Gamma.hitAvgPower = mean(cat(4,Gamma.hitPower{:}),4);
Gamma.hitAvgPhase = mean(cat(4,Gamma.hitPhase{:}),4);
Beta.missAvgPower = mean(cat(4,Beta.missPower{:}),4);
Beta.missAvgPhase = mean(cat(4,Beta.missPhase{:}),4);
Gamma.missAvgPower = mean(cat(4,Gamma.missPower{:}),4);
Gamma.missAvgPhase = mean(cat(4,Gamma.hitPhase{:}),4);

figure();
subplot(2,1,1);hold on;title('Average Beta band power for Hits');xlabel('Time (s)');ylabel('Power (in dB)');
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.RT,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
plot(IntanBehaviour.cueHitTrace(1).time,10*log10(reshape(Beta.hitAvgPower,[],size(Beta.hitAvgPower,3))),'Color', [0 0 1 0.1],'LineWidth',1.5);
h(1)=plot(IntanBehaviour.cueHitTrace(1).time,10*log10(mean(reshape(Beta.hitAvgPower,[],size(Beta.hitAvgPower,3)),1)),'Color', [0 0 1 1],'LineWidth',1.5);
plot(IntanBehaviour.cueMissTrace(1).time,10*log10(reshape(Beta.missAvgPower,[],size(Beta.missAvgPower,3))),'Color', [1 0 0 0.1],'LineWidth',1.5);
h(2)=plot(IntanBehaviour.cueMissTrace(1).time,10*log10(mean(reshape(Beta.missAvgPower,[],size(Beta.missAvgPower,3)),1)),'Color', [1 0 0 1],'LineWidth',1.5);
legend(h,'Hits','Misses');
subplot(2,1,2);hold on;title('Average Gamma band power for Hits');xlabel('Time (s)');ylabel('Power (in dB)');
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.RT,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top')
plot(IntanBehaviour.cueHitTrace(1).time,10*log10(reshape(Gamma.hitAvgPower,[],size(Gamma.hitAvgPower,3))),'Color', [0 0 1 0.1],'LineWidth',1.5);
plot(IntanBehaviour.cueHitTrace(1).time,10*log10(mean(reshape(Gamma.hitAvgPower,[],size(Gamma.hitAvgPower,3)),1)),'Color', [0 0 1 1],'LineWidth',1.5);
plot(IntanBehaviour.cueMissTrace(1).time,10*log10(reshape(Gamma.missAvgPower,[],size(Gamma.missAvgPower,3))),'Color', [1 0 0 0.1],'LineWidth',1.5);
plot(IntanBehaviour.cueMissTrace(1).time,10*log10(mean(reshape(Gamma.missAvgPower,[],size(Gamma.missAvgPower,3)),1)),'Color', [1 0 0 1],'LineWidth',1.5);

% z-scoring 
Beta.muHit = mean(cat(4,Beta.hitPower{:}),[3,4]);
Beta.sigmaHit = std(cat(4,Beta.hitPower{:}),0,[3,4]);
Beta.hitAvgPowerz = (Beta.hitAvgPower - Beta.muHit)./Beta.sigmaHit;
Beta.muMiss = mean(cat(4,Beta.missPower{:}),[3,4]);
Beta.sigmaMiss = std(cat(4,Beta.missPower{:}),0,[3,4]);
Beta.missAvgPowerz = (Beta.missAvgPower - Beta.muMiss)./Beta.sigmaMiss;
Gamma.muHit = mean(cat(4,Gamma.hitPower{:}),[3,4]);
Gamma.sigmaHit = std(cat(4,Gamma.hitPower{:}),0,[3,4]);
Gamma.hitAvgPowerz = (Gamma.hitAvgPower - Gamma.muHit)./Gamma.sigmaHit;
Gamma.muMiss = mean(cat(4,Gamma.missPower{:}),[3,4]);
Gamma.sigmaMiss = std(cat(4,Gamma.missPower{:}),0,[3,4]);
Gamma.missAvgPowerz = (Gamma.missAvgPower - Gamma.muMiss)./Gamma.sigmaMiss;

figure();
subplot(2,1,1);hold on;title('z-scored average Beta band power for Hits');xlabel('Time (s)');ylabel('Power (in dB)');
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.RT,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
plot(IntanBehaviour.cueHitTrace(1).time,reshape(Beta.hitAvgPowerz,[],size(Beta.hitAvgPowerz,3)),'Color', [0 0 1 0.1],'LineWidth',1.5);
h(1)=plot(IntanBehaviour.cueHitTrace(1).time,mean(reshape(Beta.hitAvgPowerz,[],size(Beta.hitAvgPowerz,3)),1),'Color', [0 0 1 1],'LineWidth',1.5);
plot(IntanBehaviour.cueMissTrace(1).time,reshape(Beta.missAvgPowerz,[],size(Beta.missAvgPowerz,3)),'Color', [1 0 0 0.1],'LineWidth',1.5);
h(2)=plot(IntanBehaviour.cueMissTrace(1).time,mean(reshape(Beta.missAvgPowerz,[],size(Beta.missAvgPowerz,3)),1),'Color', [1 0 0 1],'LineWidth',1.5);
legend(h,'Hits','Misses');
subplot(2,1,2);title('z-scored average Gamma band power for Hits');xlabel('Time (s)');ylabel('Power (in dB)');
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.RT,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
plot(IntanBehaviour.cueHitTrace(1).time,reshape(Gamma.hitAvgPowerz,[],size(Gamma.hitAvgPowerz,3)),'Color', [0 0 1 0.1],'LineWidth',1.5);
hold on; plot(IntanBehaviour.cueHitTrace(1).time,mean(reshape(Gamma.hitAvgPowerz,[],size(Gamma.hitAvgPowerz,3)),1),'Color', [0 0 1 1],'LineWidth',1.5);
plot(IntanBehaviour.cueMissTrace(1).time,reshape(Gamma.missAvgPowerz,[],size(Gamma.missAvgPowerz,3)),'Color', [1 0 0 0.1],'LineWidth',1.5);
plot(IntanBehaviour.cueMissTrace(1).time,mean(reshape(Gamma.missAvgPowerz,[],size(Gamma.missAvgPowerz,3)),1),'Color', [1 0 0 1],'LineWidth',1.5);


%% 
b = 16;
% Getting Amplitude of waves
wavePhases = angle(Waves.wavesHit(2).p(:,:,Waves.wavesHit(2).waveTime(b,1):Waves.wavesHit(2).waveTime(b,2)));
phi = mapAngle360(rad2deg(reshape(wavePhases,32,size(wavePhases,3))));

waveAmp = abs(Waves.wavesHit(2).p(:,:,Waves.wavesHit(2).waveTime(b,1):Waves.wavesHit(2).waveTime(b,2)));
amp = reshape(waveAmp,32,size(waveAmp,3));

figure;
subplot(2,1,1);
plot(phi');
subplot(2,1,2);
plot(amp');

xijtHit = cellfun(@(s) s(1,1,1),xgpHit);
xijtMiss = cellfun(@(s) s(1,1,1),xgpMiss);

infomap = reshape(MIPhase,[],size(MIPhase,3));
infopeak = mean(infomap,2);
[~,indexsort] = sort(infopeak,1,'descend');

sortedInfomap = infomap(indexsort,:);

figure();
title("Mututal Information across all electrodes - Phase")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,sortedInfomap); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Information (bits)';
xline(0,'-w','Cue','LabelVerticalAlignment','top');

figure();
for i=1:size(IntanBehaviour.missTrace,2)
    plot(IntanBehaviour.missTrace(i).time,IntanBehaviour.missTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.missTrace(1).time,mean(horzcat(IntanBehaviour.missTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Misses');box off;

%% Performing wave stats for a time window

%% Phase allignement but for phase gradient direction
xgp = arrayfun(@(s) s.pd, Waves.wavesHit, 'UniformOutput', false);
PGAHit = getPhaseAlignment(xgp,parameters);
xgp = arrayfun(@(s) s.pd, Waves.wavesMiss, 'UniformOutput', false);
PGAMiss = getPhaseAlignment(xgp,parameters);

figure();
subplot(2,1,1);
title("Phase Alignment across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PGAHit,[],size(PGAHit,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Phase Alignment across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PGAMiss,[],size(PGAMiss,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');


figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PGAHit,[1 2])),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PGAMiss,[1 2])),'-k','LineWidth',1);
ylabel("Phase Gradient Alignment"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.RT,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Phase Gradient Alignment for Hits');box off;% legend('Hits','Miss');

%%
figure();
for i = 1:4
    st = (i-1)*500;
    sp = (i)*500;
    dirComb1 = horzcat(Waves.wavesHit(1:end).waveDir);
    evalPointsHit = horzcat(Waves.wavesHit(1:end).evaluationPoints);
    WaveComb(i).Hit = dirComb1(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombMiss1 = horzcat(Waves.wavesMiss(1:end).waveDir);
    evalPointsMiss = horzcat(Waves.wavesMiss(1:end).evaluationPoints);
    WaveComb(i).Miss = dirCombMiss1(evalPointsMiss >=st & evalPointsMiss <= sp);

    subplot(2,4,i)
    polarhistogram(WaveComb(i).Hit,60);
    subplot(2,4,4+i)
    polarhistogram(WaveComb(i).Miss,60);
end

figure(); hold on;
plot(1:4,arrayfun(@(s) mean(s.Hit,'all'), WaveComb));
plot(1:4,arrayfun(@(s) mean(s.Miss,'all'), WaveComb));
a = arrayfun(@(s) circ_kuipertest(s.Hit,s.Miss,60,0) , WaveComb);
%%

for i = 1:30
    st = (i-1)*100;
    sp = (i)*100;
    dirComb1 = horzcat(Waves.wavesHit(1:end).speed);
    evalPointsHit = horzcat(Waves.wavesHit(1:end).evaluationPoints);
    WaveComb(i).Hit = dirComb1(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombMiss1 = horzcat(Waves.wavesMiss(1:end).speed);
    evalPointsMiss = horzcat(Waves.wavesMiss(1:end).evaluationPoints);
    WaveComb(i).Miss = dirCombMiss1(evalPointsMiss >=st & evalPointsMiss <= sp);
end

figure(); hold on;
avgVarHit = arrayfun(@(s) mean(s.Hit,'all'), WaveComb);
maxVarHit = max(avgVarHit);
avgVarMiss = arrayfun(@(s) mean(s.Miss,'all'), WaveComb);
maxVarMiss = max(avgVarMiss);
a = arrayfun(@(s) ranksum(s.Hit,s.Miss) , WaveComb);
pVal = (a<0.01);
significanceX = find(pVal)*100;
significanceY = 1.1*max(maxVarHit,maxVarMiss)*ones(1,length(significanceX));
h1 = plot(100:100:3000,(avgVarHit),'Color', [0 0 0 1],'LineWidth',1.5);
h2 = plot(100:100:3000,(avgVarMiss),'Color', [0 0 0 0.2],'LineWidth',1.5);
xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Wave speed (cm/s)');
plot(significanceX,significanceY,'r*');
legend([h1 h2],'Hits','Miss','Location','best');
xlim([1000 3000]);set(gca,'TickDir','out','fontsize',14')

% figure,plot(100:100:3000,a);

%%
for i = 1:30
    st = (i-1)*100;
    sp = (i)*100;
    dirComb1 = horzcat(Waves.wavesHitReward(1:end).speed);
    evalPointsHit = horzcat(Waves.wavesHitReward(1:end).evaluationPoints);
    WaveComb(i).Hit = dirComb1(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombMiss1 = horzcat(Waves.wavesFA(1:end).speed);
    evalPointsMiss = horzcat(Waves.wavesFA(1:end).evaluationPoints);
    WaveComb(i).Miss = dirCombMiss1(evalPointsMiss >=st & evalPointsMiss <= sp);
end

figure(); hold on;
avgVarHit = arrayfun(@(s) mean(s.Hit,'all'), WaveComb);
maxVarHit = max(avgVarHit);
avgVarMiss = arrayfun(@(s) mean(s.Miss,'all'), WaveComb);
maxVarMiss = max(avgVarMiss);
a = arrayfun(@(s) ranksum(s.Hit,s.Miss) , WaveComb);
pVal = (a<0.05);
significanceX = find(pVal)*100;
significanceY = 1.1*max(maxVarHit,maxVarMiss)*ones(1,length(significanceX));
h1 = plot(100:100:3000,avgVarHit,'Color', [0 0 0 1],'LineWidth',1.5);
h2 = plot(100:100:3000,avgVarMiss,'Color', [0 0 0 0.2],'LineWidth',1.5);
xline(1501,'--r','Reward','LabelVerticalAlignment','top');
ylabel('Speed (cm/s)');xlabel('Time (ms)')
plot(significanceX,significanceY,'r*');
legend([h1 h2],'Hits','FAs');

% figure,plot(100:100:3000,a);

%%
for i = 1:30
    st = (i-1)*100;
    sp = (i)*100;
    dirComb1 = horzcat(Waves.wavesMIHit(1:end).speed);
    evalPointsHit = horzcat(Waves.wavesMIHit(1:end).evaluationPoints);
    WaveComb(i).Hit = dirComb1(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombMiss1 = horzcat(Waves.wavesMIFA(1:end).speed);
    evalPointsMiss = horzcat(Waves.wavesMIFA(1:end).evaluationPoints);
    WaveComb(i).Miss = dirCombMiss1(evalPointsMiss >=st & evalPointsMiss <= sp);
end

figure(); hold on;
avgVarHit = arrayfun(@(s) mean(s.Hit,'all'), WaveComb);
maxVarHit = max(avgVarHit);
avgVarMiss = arrayfun(@(s) mean(s.Miss,'all'), WaveComb);
maxVarMiss = max(avgVarMiss);
a = arrayfun(@(s) ranksum(s.Hit,s.Miss) , WaveComb);
pVal = (a<0.05);
significanceX = find(pVal)*100;
significanceY = 1.1*max(maxVarHit,maxVarMiss)*ones(1,length(significanceX));
h1 = plot(100:100:3000,avgVarHit,'Color', [0 0 0 1],'LineWidth',1.5);
h2 = plot(100:100:3000,avgVarMiss,'Color', [0 0 0 0.2],'LineWidth',1.5);
xline(1501,'--r','MI');xlabel('Time (ms)')
plot(significanceX,significanceY,'r*');
legend([h1 h2],'Hits','FAs');
set(gca,'TickDir','out','fontsize',14')
xlim([1000 3000]);set(gca,'TickDir','out','fontsize',14');ylabel('Wave speed (cm/s)');
% figure,plot(100:100:3000,a);


%%
for i = 1:6
    st = (i-1)*500;
    sp = (i)*500;
    dirComb1 = horzcat(Waves.wavesHit(1:end).waveDir);
    evalPointsHit = horzcat(Waves.wavesHit(1:end).evaluationPoints);
    WaveCombDir(i).Hit = dirComb1(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombMiss1 = horzcat(Waves.wavesMiss(1:end).waveDir);
    evalPointsMiss = horzcat(Waves.wavesMiss(1:end).evaluationPoints);
    WaveCombDir(i).Miss = dirCombMiss1(evalPointsMiss >=st & evalPointsMiss <= sp);
end

figure(); hold on;
for i=1:6
    subplot(2,6,i);
    polarhistogram(WaveCombDir(i).Hit,60);
    subplot(2,6,6+i);
    polarhistogram(WaveCombDir(i).Miss,60);
end

circpVal = arrayfun(@(s)  circ_kuipertest(s.Hit,s.Miss,100,0) , WaveCombDir);

for i = 1:6
    st = (i-1)*500;
    sp = (i)*500;
    dirComb1 = horzcat(Waves.wavesHitReward(1:end).waveDir);
    evalPointsHit = horzcat(Waves.wavesHitReward(1:end).evaluationPoints);
    WaveCombDir(i).Hit = dirComb1(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombMiss1 = horzcat(Waves.wavesFA(1:end).waveDir);
    evalPointsMiss = horzcat(Waves.wavesFA(1:end).evaluationPoints);
    WaveCombDir(i).Miss = dirCombMiss1(evalPointsMiss >=st & evalPointsMiss <= sp);
end

figure(); hold on;
for i=1:6
    subplot(2,6,i);
    polarhistogram(WaveCombDir(i).Hit,60);
    subplot(2,6,6+i);
    polarhistogram(WaveCombDir(i).Miss,60);
end
circpVal = arrayfun(@(s)  circ_kuipertest(s.Hit,s.Miss,100,0) , WaveCombDir);


for i = 1:6
    st = (i-1)*500;
    sp = (i)*500;
    dirComb1 = horzcat(Waves.wavesMIHit(1:end).waveDir);
    evalPointsHit = horzcat(Waves.wavesMIHit(1:end).evaluationPoints);
    WaveCombDir(i).Hit = dirComb1(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombMiss1 = horzcat(Waves.wavesMIFA(1:end).waveDir);
    evalPointsMiss = horzcat(Waves.wavesMIFA(1:end).evaluationPoints);
    WaveCombDir(i).Miss = dirCombMiss1(evalPointsMiss >=st & evalPointsMiss <= sp);
end

figure(); hold on;
for i=1:6
    subplot(2,6,i);
    polarhistogram(WaveCombDir(i).Hit,60);
    subplot(2,6,6+i);
    polarhistogram(WaveCombDir(i).Miss,60);
end
circpVal = arrayfun(@(s)  circ_kuipertest(s.Hit,s.Miss,100,0) , WaveCombDir);

%% test for direction
dirComb = horzcat(Waves.wavesHit(1:end).waveDir);
dirCombMiss = horzcat(Waves.wavesMiss(1:end).waveDir);

[pval, ~, ~] = circ_kuipertest(dirComb, dirCombMiss, 100, 0);

%%
a = horzcat(Waves.wavesMIHit(:).nWaves);
b = horzcat(Waves.wavesMIFA(:).nWaves);

p = ranksum(a,b)

figure('Name','Number of detected waves   Hits and Misses');
group = [ones(size(a')); 2.*ones(size(b'))];
boxplot([a';b'],group,'BoxStyle','filled','PlotStyle','compact');box off;
set(gca,'XTickLabel',{'Hits','Misses'});
xlabel('Number of waves');


%%
figure;
hold on;
arrayfun(@(s) plot(s.time, s.trace),IntanBehaviour.cueMissTrace,'UniformOutput',false);

%% Plotting waves 
a = Waves.wavesHit(1).xf(:,:,1100:1200);
a = reshape(a,[],size(a,3));
[peak, peakIndex] = max(a,[],2);

[~, sortindx] = sort(peakIndex);
asorted = a(sortindx,:);

figure();
plot(asorted','color', [.5 .5 .5], 'linewidth', 1.5 );
xlabel( 'Time (ms)' ); ylabel( 'Amplitude (\muV)' );


plot_evaluation_points( Waves.wavesHit(14).p, Waves.wavesHit(14).evaluationPoints );


%% 
animateWaves(42,IntanBehaviour.cueHitTrace,Waves.wavesHit,0,13);

figure,stack_plot(reshape(IntanBehaviour.cueHitTrace(23).xf,[],3001),0,4,1000);

figure,stack_plot(reshape(IntanBehaviour.hitTrace(23).xf,[],3001),1,4,1000);

figure,stack_plot(reshape(LFP.xf(:,:,1:10000),[],10000),0,4,1000);


%%
wavesHitPresent = vertcat(Waves.wavesHit.wavePresent);
wavesMissPresent = vertcat(Waves.wavesMiss.wavePresent);
wavesHitStart = vertcat(Waves.wavesHit.waveStart);
wavesMissStart = vertcat(Waves.wavesMiss.waveStart);

RTTraceTime = (IntanBehaviour.reactionTime*parameters.Fs) +(parameters.windowBeforeCue*parameters.Fs);

figure();
subplot(4,1,1);
title('Waves During Hit Trials')
rasterPlot(wavesHitPresent);hold on;
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
rasterPlot(wavesMissPresent);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);%ylim([1 20]);
set(gca,'TickDir','out','fontsize',14'); box off;
subplot(4,1,4)
bar((sum(wavesMissPresent,1)/size(IntanBehaviour.cueMissTrace,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Cue','LabelVerticalAlignment','top');
ylim([0 0.3]);xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs])
ylabel('Wave probability');xlabel('Time (in ms)')
set(gca,'TickDir','out','fontsize',14'); box off;

%%
wavesHitRewardPresent = vertcat(Waves.wavesHitReward.wavePresent);
wavesFAPresent = vertcat(Waves.wavesFA.wavePresent);
wavesHitRewardStart = vertcat(Waves.wavesHitReward.waveStart);
wavesFAStart = vertcat(Waves.wavesFA.waveStart);

figure();
subplot(4,1,1);
title('Waves During Hit Trials')
rasterPlot(wavesHitRewardPresent);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Reward','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([1 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);
subplot(4,1,2)
bar((sum(wavesHitRewardPresent,1)/size(IntanBehaviour.hitTrace,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','Reward','LabelVerticalAlignment','top');
ylim([0 0.2]);
ylabel('Wave probability');xlabel('Time (in ms)')
subplot(4,1,3);
title('Waves During FA Trials')
rasterPlot(wavesFAPresent);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','No Reward','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([1 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);
subplot(4,1,4)
bar((sum(wavesFAPresent,1)/size(IntanBehaviour.missTrace,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','No Reward','LabelVerticalAlignment','top');
ylim([0 0.2]);
ylabel('Wave probability');xlabel('Time (in ms)')

%%

wavesMIHitPresent = vertcat(Waves.wavesMIHit.wavePresent);
wavesMIFAPresent = vertcat(Waves.wavesMIFA.wavePresent);
wavesMIHitStart = vertcat(Waves.wavesMIHit.waveStart);
wavesMIFAStart = vertcat(Waves.wavesMIFA.waveStart);

figure();
subplot(4,1,1);
title('Waves During Hit ')
rasterPlot(wavesMIHitPresent);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
xline(-(mean(IntanBehaviour.reactionTime,'all')*parameters.Fs + parameters.windowBeforeCue*parameters.Fs+1),'--m','Avg. Cue','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);ylim([1 80]);
set(gca,'TickDir','out','fontsize',14'); box off;
subplot(4,1,2)
bar((sum(wavesMIHitPresent,1)/size(IntanBehaviour.MIHitTrace,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
xline(-(mean(IntanBehaviour.reactionTime,'all')*parameters.Fs + parameters.windowBeforeCue*parameters.Fs+1),'--m','Avg. Cue','LabelVerticalAlignment','top');
ylim([0 0.3]);xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);
ylabel('Wave probability');xlabel('Time (in ms)');set(gca,'TickDir','out','fontsize',14'); box off;
subplot(4,1,3);
title('Waves During FA')
rasterPlot(wavesMIFAPresent);hold on;
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
ylabel('Trials');xlabel('Time (in ms)');xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs]);ylim([1 80]);
set(gca,'TickDir','out','fontsize',14'); box off;
subplot(4,1,4)
bar((sum(wavesMIFAPresent,1)/size(IntanBehaviour.MIFATrace,2)));
xline(parameters.windowBeforeCue*parameters.Fs+1,'--r','MI','LabelVerticalAlignment','top');
ylim([0 0.3]);xlim([1000 parameters.windowBeforeCue*parameters.Fs+1+parameters.windowAfterCue*parameters.Fs])
ylabel('Wave probability');xlabel('Time (in ms)')
set(gca,'TickDir','out','fontsize',14'); box off;


%%
LFP.xf(4,2,:) = NaN;
LFP.xgp(4,2,:) = NaN;
LFP.wt(4,2,:) = NaN;
LFP.xf(5,1,:) = NaN;
LFP.xgp(5,1,:) = NaN;
LFP.wt(5,1,:) = NaN;

%% 
m = matfile('D:\Om\GridsShanksDay4.mat');
fpath = m.fpath;
parameters = m.parameters;

IntanBehaviour.cueHitTrace = m.IntanBehaviour.cueHitTrace;

%%
betaBurstAllHit = vertcat(BetaEventHit.betaBurstPresent);
betaBurstAllMiss = vertcat(BetaEventMiss.betaBurstPresent);
figure();
subplot(411);
rasterPlot(betaBurstAllHit);
xline(1501,'--r','MI','LabelVerticalAlignment','top');
% xline((mean(IntanBehaviour.RT,'all')*parameters.Fs + 1501),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Beta Burst');xlabel('Time (ms)');xlim([1 3001]);
subplot(412);
bar(sum(betaBurstAllHit,1)/size(IntanBehaviour.cueHitTrace,2));
ylabel('Beta Burst Probability');xlabel('Time (ms)');
xline(1501,'--r','MI','LabelVerticalAlignment','top');
subplot(413);
rasterPlot(betaBurstAllMiss);
xline(1501,'--r','MI','LabelVerticalAlignment','top');
ylabel('Beta Burst');xlabel('Time (ms)');xlim([1 3001]);
subplot(414);
bar(sum(betaBurstAllMiss,1)/size(IntanBehaviour.cueMissTrace,2));
ylabel('Beta Burst Probability');xlabel('Time (ms)');
xline(1501,'--r','MI','LabelVerticalAlignment','top');


betaBurstAllHit = vertcat(GammaEventHit.betaBurstPresent);
betaBurstAllMiss = vertcat(GammaEventMiss.betaBurstPresent);

figure();
subplot(411);
rasterPlot(betaBurstAllHit);
xline(1501,'--r','MI','LabelVerticalAlignment','top');
% xline((mean(IntanBehaviour.reactionTime,'all')*parameters.Fs + 1501),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Gamma Burst');xlabel('Time (ms)');xlim([1 3001]);
subplot(412);
bar(sum(betaBurstAllHit,1)/size(IntanBehaviour.MIHitTrace,2));
ylabel('Gamma Burst Probability');xlabel('Time (ms)');
xline(1501,'--r','MI','LabelVerticalAlignment','top');
subplot(413);
rasterPlot(betaBurstAllMiss);
xline(1501,'--r','MI','LabelVerticalAlignment','top');
ylabel('Gamma Burst');xlabel('Time (ms)');xlim([1 3001]);
subplot(414);
bar(sum(betaBurstAllMiss,1)/size(IntanBehaviour.MIFATrace,2));
ylabel('Gamma Burst Probability');xlabel('Time (ms)');
xline(1501,'--r','MI','LabelVerticalAlignment','top');

%% Granger Causality
datGrid = vertcat(cell2mat(arrayfun(@(x) reshape(x.rawLFP,32,[]), IntanBehaviour.cueHitTrace, 'UniformOutput', false)));
datGrid = removeNaNRows(datGrid);
datGrid = reshape(datGrid,30,3001,[]);
datProbe = vertcat(cell2mat(arrayfun(@(x) reshape(x.rawLFPProbe,64,[]), IntanBehaviour.cueHitTrace, 'UniformOutput', false)));
datProbe = removeNaNRows(datProbe);
datProbe = reshape(datProbe,64,3001,[]);

dat =vertcat(datGrid,datProbe);
GCM22M1 = mvgcstate(dat); %dat structured as (node,samples,trials)
%%
 % Plot spectral causal graph.
 figure(); clf;
 sgtitlex('Pairwise-conditional Granger causality - frequency domain');
 plot_spw(output.P(1:3,1:3,:),output.ops.fs);


%%
for i = 1:20
    st = (i-1)*100;
    sp = (i)*100;
    dirComb1 = horzcat(Waves.wavesHitReward(1:end).speed);
    evalPointsHit = horzcat(Waves.wavesHitReward(1:end).evaluationPoints);
    WaveComb(i).Hit = dirComb1(evalPointsHit >=st & evalPointsHit <= sp);

    dirCombMiss1 = horzcat(Waves.wavesFA(1:end).speed);
    evalPointsMiss = horzcat(Waves.wavesFA(1:end).evaluationPoints);
    WaveComb(i).Miss = dirCombMiss1(evalPointsMiss >=st & evalPointsMiss <= sp);
end

% figure(); 
subplot(1,2,1)
hold on;
plot(100:100:2000,arrayfun(@(s) mean(s.Hit,'all'), WaveComb),'Color', [0 0 0 1],'LineWidth',1.5);
plot(100:100:2000,arrayfun(@(s) mean(s.Miss,'all'), WaveComb),'Color', [0 0 0 0.2],'LineWidth',1.5);
xline(501,'--r','Cue','LabelVerticalAlignment','top');
legend('Hits','Misses');
ylabel('Speed (cm/s)');xlabel('Time (ms)')


a = arrayfun(@(s) ranksum(s.Hit,s.Miss) , WaveComb);
figure,plot(100:100:2000,a);

%%
IntanBehaviour.reactionTime = arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace);

figure,histogram(IntanBehaviour.reactionTime,20);

weirdRTIndex = find(IntanBehaviour.reactionTime<0 | IntanBehaviour.reactionTime > 2);

IntanBehaviour.reactionTime(weirdRTIndex) = [];
IntanBehaviour.cueHitTrace(weirdRTIndex) = [];
IntanBehaviour.hitTrace(weirdRTIndex) = [];
IntanBehaviour.MIHitTrace(weirdRTIndex) = [];
Waves.wavesHit(weirdRTIndex) = [];
Waves.wavesMIHit(weirdRTIndex) = [];
Waves.wavesHitReward(weirdRTIndex) = [];


%% 

xgp = reshape(IntanBehaviour.cueHitTrace(60).xgp,32,[]);
phase = angle(xgp);
figure,stack_plot(phase,0,10,1000);hold on; xline(1.501)

PR = zeros(size(IntanBehaviour.cueHitTrace,2),size(IntanBehaviour.cueHitTrace(1).xgp,3));
for i=1:size(IntanBehaviour.cueHitTrace,2)
    xgp = reshape(IntanBehaviour.cueHitTrace(i).xgp,32,[]);
    xgpnorm = xgp./abs(xgp);
    xgpnorm = removeNaNRows(xgpnorm);     
    for j=1:size(xgpnorm,2)
        PR(i,j) = abs(sum(xgpnorm(:,j))/size(xgpnorm,1));
    end
end

%% 
avgSpeedRT = zeros(1,size(IntanBehaviour.cueHitTrace,2));
for i=1:size(IntanBehaviour.cueHitTrace,2)
    evalpoints = Waves.wavesHit(i).evaluationPoints;
    a = find(evalpoints>1500 & evalpoints<1700);
    % b = find(evalpoints>0000 & evalpoints<1500);
    avgSpeedRT(i) = mean(Waves.wavesHit(i).speed(a),'all','omitnan');% - mean(Waves.wavesHit(i).speed(b),'all','omitnan');;
end

mdl = fitlm(avgSpeedRT,IntanBehaviour.reactionTime)
figure,plot(mdl);

%%
firstWave = zeros(1,size(IntanBehaviour.cueHitTrace,2));
for i=1:size(IntanBehaviour.cueHitTrace,2)
    a = find(Waves.wavesHit(i).waveStart(1501:1701),1);
    if isempty(a)
        firstWave(i) = NaN;
    else
        firstWave(i) = a;
    end
end

mdl = fitlm(firstWave,IntanBehaviour.reactionTime)
figure,plot(mdl);


%% 
pgdTrial = vertcat(Waves.wavesHit.PGD);
figure,imagesc(IntanBehaviour.cueHitTrace(1).time,1:size(Waves.wavesHit,2),pgdTrial);
colormap(hot);axis xy;


%%
i = 1;
j = 1;
badIndex = [];
while i<=121
    if IntanBehaviour.cueHitTrace(i).LFPIndex(1) == cueHitTraceLFP(j).LFPIndex(1)
        i = i+1;
        j = j+1;
    else
        badIndex = [badIndex;j];
        j = j+1;
    end
end

cueHitTraceLFP(badIndex) = [];


for i=1:121
    IntanBehaviour.cueHitTrace(i).rawLFPProbe = cueHitTraceLFP(i).rawLFPProbe;
    IntanBehaviour.cueHitTrace(i).xfProbe = cueHitTraceLFP(i).xfProbe;
    IntanBehaviour.cueHitTrace(i).xgpProbe = cueHitTraceLFP(i).xgpProbe;
    IntanBehaviour.cueHitTrace(i).wtProbe = cueHitTraceLFP(i).wtProbe;
end


for i=1:70
    IntanBehaviour.cueMissTrace(i).rawLFPProbe = cueMissTraceLFP(i).rawLFPProbe;
    IntanBehaviour.cueMissTrace(i).xfProbe = cueMissTraceLFP(i).xfProbe;
    IntanBehaviour.cueMissTrace(i).xgpProbe = cueMissTraceLFP(i).xgpProbe;
    IntanBehaviour.cueMissTrace(i).wtProbe = cueMissTraceLFP(i).wtProbe;
end

%%
i = 1;
j = 1;
badIndex = [];
while i<=120
    if IntanBehaviour.MIHitTrace(i).LFPIndex(1) == MIHitTraceProbe(j).LFPIndex(1)
        i = i+1;
        j = j+1;
    else
        badIndex = [badIndex;j];
        j = j+1;
    end
end

MIHitTraceLFP(badIndex) = [];


for i=1:121
    IntanBehaviour.MIHitTrace(i).rawLFPProbe = MIHitTraceLFP(i).rawLFPProbe;
    IntanBehaviour.MIHitTrace(i).xfProbe = MIHitTraceLFP(i).xfProbe;
    IntanBehaviour.MIHitTrace(i).xgpProbe = MIHitTraceLFP(i).xgpProbe;
    IntanBehaviour.MIHitTrace(i).wtProbe = MIHitTraceLFP(i).wtProbe;
end

%%
hitProbeXGP = cell(1,size(IntanBehaviour.cueHitTrace,2));
for i=1:size(IntanBehaviour.cueHitTrace,2)
    probeAngle = angle(IntanBehaviour.cueHitTrace(i).xgpProbe(linearProbe,:));
    gridAngle = angle(squeeze(IntanBehaviour.cueHitTrace(i).xgp(5,3,:)))';
    hitProbeXGP{1,i} = reshape(exp(1i*(probeAngle-repmat(gridAngle,length(linearProbe),1))),length(linearProbe),1,[]);
end

missProbeXGP = cell(1,size(IntanBehaviour.cueMissTrace,2));
for i=1:size(IntanBehaviour.cueMissTrace,2)
    probeAngle = angle(IntanBehaviour.cueMissTrace(i).xgpProbe(linearProbe,:));
    gridAngle = angle(squeeze(IntanBehaviour.cueMissTrace(i).xgp(5,3,:)))';
    missProbeXGP{1,i} = reshape(exp(1i*(probeAngle-repmat(gridAngle,length(linearProbe),1))),length(linearProbe),1,[]);
end
MIHitProbeXGP = cell(1,size(IntanBehaviour.MIHitTrace,2));
for i=1:size(IntanBehaviour.MIHitTrace,2)
    probeAngle = angle(IntanBehaviour.MIHitTrace(i).xgpProbe(linearProbe,:));
    gridAngle = angle(squeeze(IntanBehaviour.MIHitTrace(i).xgp(5,2,:)))';
    MIHitProbeXGP{1,i} = reshape(exp(1i*(probeAngle-repmat(gridAngle,length(linearProbe),1))),length(linearProbe),1,[]);
end

MIFAProbeXGP = cell(1,size(IntanBehaviour.MIFATrace,2));
for i=1:size(IntanBehaviour.MIFATrace,2)
    probeAngle = angle(IntanBehaviour.MIFATrace(i).xgpProbe(linearProbe,:));
    gridAngle = angle(squeeze(IntanBehaviour.MIFATrace(i).xgp(5,2,:)))';
    MIFAProbeXGP{1,i} = reshape(exp(1i*(probeAngle-repmat(gridAngle,length(linearProbe),1))),length(linearProbe),1,[]);
end

ISPC.Hit = calPhaseAlignment(hitProbeXGP,parameters);
ISPC.Miss = calPhaseAlignment(missProbeXGP,parameters);
ISPC.Hit = calPhaseAlignment(MIHitProbeXGP,parameters);
ISPC.Miss = calPhaseAlignment(miProbeXGP,parameters);

figure();
title("Inter-site Phase Clustering")
subplot(2,1,1);
imagesc(IntanBehaviour.cueHitTrace(1).time,1:length(linearProbe),reshape(ISPC.Hit,[],size(ISPC.Hit,3))); colormap(sky);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');clim([0.45 1]);
subplot(2,1,2);
title("Phase Alignment across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:length(linearProbe),reshape(ISPC.Miss,[],size(ISPC.Miss,3))); colormap(sky);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');clim([0.45 1]);


%% 

IntanBehaviour1 = IntanBehaviour;
IntanBehaviour1 = rmfield(IntanBehaviour1,{'leverTrace','time','rewardTrace','cueTrace','nCueHit','nCueMiss','cueHit','cueMiss','nHit','nMiss'});
IntanBehaviour1.cueHitTrace = rmfield(IntanBehaviour1.cueHitTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.cueMissTrace = rmfield(IntanBehaviour1.cueMissTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.hitTrace = rmfield(IntanBehaviour1.hitTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.missTrace = rmfield(IntanBehaviour1.missTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.MIHitTrace = rmfield(IntanBehaviour1.MIHitTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.MIFATrace = rmfield(IntanBehaviour1.MIFATrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xfgamma', 'xgpgamma', 'wtgamma'});

IntanBehaviour1.cueHitTrace = rmfield(IntanBehaviour1.cueHitTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xftheta', 'xgptheta', 'wttheta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.cueMissTrace = rmfield(IntanBehaviour1.cueMissTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xftheta', 'xgptheta', 'wttheta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.hitTrace = rmfield(IntanBehaviour1.hitTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xftheta', 'xgptheta', 'wttheta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.missTrace = rmfield(IntanBehaviour1.missTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xftheta', 'xgptheta', 'wttheta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.MIHitTrace = rmfield(IntanBehaviour1.MIHitTrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xftheta', 'xgptheta', 'wttheta', 'xfgamma', 'xgpgamma', 'wtgamma'});
IntanBehaviour1.MIFATrace = rmfield(IntanBehaviour1.MIFATrace, {'rawLFP', 'xfbeta', 'xgpbeta', 'wtbeta', 'xftheta', 'xgptheta', 'wttheta', 'xfgamma', 'xgpgamma', 'wtgamma'});

