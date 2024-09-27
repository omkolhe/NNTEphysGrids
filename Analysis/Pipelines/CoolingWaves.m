%% Setting saving path 
savepath = uigetdir(fpath);
mkdir([savepath '\Figures'])

%% Spliting IntanBehaviour and waves into baseline, cooling and recovery 
baselineIndex = [];
coolingIndex = [];
recoveryIndex = [];
recoveryflag = 0;
for i = 1:size(IntanBehaviour.cueHitTrace,2)
    if IntanBehaviour.cueHitTrace(i).temp <= 12
        coolingIndex = i;
        recoveryflag =1;
    elseif IntanBehaviour.cueHitTrace(i).temp >= 20 && recoveryflag == 0
       baselineIndex = i;
    elseif IntanBehaviour.cueHitTrace(i).temp >= 15 && recoveryflag == 1
        recoveryIndex = i;
    end
end
IntanBehaviourBaseline.cueHitTrace = IntanBehaviour.cueHitTrace(1:baselineIndex);
IntanBehaviourCooling.cueHitTrace = IntanBehaviour.cueHitTrace(baselineIndex+1:coolingIndex);
IntanBehaviourRecovery.cueHitTrace = IntanBehaviour.cueHitTrace(coolingIndex+1:end);
IntanBehaviourBaseline.hitTrace = IntanBehaviour.hitTrace(1:baselineIndex);
IntanBehaviourCooling.hitTrace = IntanBehaviour.hitTrace(baselineIndex+1:coolingIndex);
IntanBehaviourRecovery.hitTrace = IntanBehaviour.hitTrace(coolingIndex+1:end);
IntanBehaviourBaseline.MIHitTrace = IntanBehaviour.MIHitTrace(1:baselineIndex);
IntanBehaviourCooling.MIHitTrace = IntanBehaviour.MIHitTrace(baselineIndex+1:coolingIndex);
IntanBehaviourRecovery.MIHitTrace = IntanBehaviour.MIHitTrace(coolingIndex+1:end);
WavesBaseline.wavesHit = Waves.wavesHit(1:baselineIndex);
WavesCooling.wavesHit = Waves.wavesHit(baselineIndex+1:coolingIndex);
WavesRecovery.wavesHit = Waves.wavesHit(coolingIndex+1:end);
WavesBaseline.wavesMIHit = Waves.wavesMIHit(1:baselineIndex);
WavesCooling.wavesMIHit = Waves.wavesMIHit(baselineIndex+1:coolingIndex);
WavesRecovery.wavesMIHit = Waves.wavesMIHit(coolingIndex+1:end);
WavesBaseline.wavesHitReward = Waves.wavesHitReward(1:baselineIndex);
WavesCooling.wavesHitReward = Waves.wavesHitReward(baselineIndex+1:coolingIndex);
WavesRecovery.wavesHitReward = Waves.wavesHitReward(coolingIndex+1:end);

baselineIndex = [];
coolingIndex = [];
recoveryIndex = [];
recoveryflag = 0;
for i = 1:size(IntanBehaviour.missTrace,2)
    if IntanBehaviour.missTrace(i).temp <= 12
        coolingIndex = i;
        recoveryflag =1;
    elseif IntanBehaviour.missTrace(i).temp >= 20 && recoveryflag == 0
       baselineIndex = i;
    elseif IntanBehaviour.missTrace(i).temp >= 15 && recoveryflag == 1
        recoveryIndex = i;
    end
end
IntanBehaviourBaseline.missTrace = IntanBehaviour.missTrace(1:baselineIndex);
IntanBehaviourCooling.missTrace = IntanBehaviour.missTrace(baselineIndex+1:coolingIndex);
IntanBehaviourRecovery.missTrace = IntanBehaviour.missTrace(coolingIndex+1:end);
IntanBehaviourBaseline.MIFATrace = IntanBehaviour.MIFATrace(1:baselineIndex);
IntanBehaviourCooling.MIFATrace = IntanBehaviour.MIFATrace(baselineIndex+1:coolingIndex);
IntanBehaviourRecovery.MIFATrace = IntanBehaviour.MIFATrace(coolingIndex+1:end);
WavesBaseline.wavesFA = Waves.wavesFA(1:baselineIndex);
WavesCooling.wavesFA = Waves.wavesFA(baselineIndex+1:coolingIndex);
WavesRecovery.wavesFA = Waves.wavesFA(coolingIndex+1:end);
WavesBaseline.wavesMIFA = Waves.wavesMIFA(1:baselineIndex);
WavesCooling.wavesMIFA = Waves.wavesMIFA(baselineIndex+1:coolingIndex);
WavesRecovery.wavesMIFA = Waves.wavesMIFA(coolingIndex+1:end);

baselineIndex = [];
coolingIndex = [];
recoveryIndex = [];
recoveryflag = 0;
for i = 1:size(IntanBehaviour.cueMissTrace,2)
    if IntanBehaviour.cueMissTrace(i).temp <= 12
        coolingIndex = i;
        recoveryflag =1;
    elseif IntanBehaviour.cueMissTrace(i).temp >= 20 && recoveryflag == 0
       baselineIndex = i;
    elseif IntanBehaviour.cueMissTrace(i).temp >= 15 && recoveryflag == 1
        recoveryIndex = i;
    end
end
IntanBehaviourBaseline.cueMissTrace = IntanBehaviour.cueMissTrace(1:baselineIndex);
IntanBehaviourCooling.cueMissTrace = IntanBehaviour.cueMissTrace(baselineIndex+1:coolingIndex);
IntanBehaviourRecovery.cueMissTrace = IntanBehaviour.cueMissTrace(coolingIndex+1:end);
WavesBaseline.wavesMiss = Waves.wavesMiss(1:baselineIndex);
WavesCooling.wavesMiss = Waves.wavesMiss(baselineIndex+1:coolingIndex);
WavesRecovery.wavesMiss = Waves.wavesMiss(coolingIndex+1:end);

%% Plotting Behaviour 

%Plotting hits vs cooling 
hitTime = cell2mat(arrayfun(@(s) s.LFPtime(1501), IntanBehaviour.cueHitTrace, 'UniformOutput', false));
RT = cell2mat(arrayfun(@(s) s.reactionTime, IntanBehaviour.cueHitTrace, 'UniformOutput', false));
h = figure();
plot(IntanBehaviour.time/60,lowpass(IntanBehaviour.tempTrace,0.1,parameters.Fs));
ylabel('Temperature (in $^\circ$ C)','Interpreter','latex');
xlabel('Time (in min)')
hold on; yyaxis right; box off;
plot(hitTime/60,RT*1000,'r*');
ylabel('Reaction Time (in ms)');
saveas(h,[savepath '\Figures\RTvsCooling.png']);
saveas(h,[savepath '\Figures\RTvsCooling.fig']);
%% Plotting reaction time 
baselineRT = [];
coolRT = [];
recoveryRT = [];
recoveryFlag = 0; % Run before loop; recoveryFlag = 0 - Baseline, = 1 - Recovery
% only works for 1 baseline, cool and recovery cycle
for i=1:size(IntanBehaviour.cueHitTrace,2)  
    if IntanBehaviour.cueHitTrace(i).temp <= 15
        coolRT = [coolRT;IntanBehaviour.cueHitTrace(i).reactionTime];
        recoveryFlag = 1;
%         recoveryFlag = 0;
    elseif IntanBehaviour.cueHitTrace(i).temp >= 20 && recoveryFlag == 0
        baselineRT= [baselineRT;IntanBehaviour.cueHitTrace(i).reactionTime];
    elseif IntanBehaviour.cueHitTrace(i).temp >= 18 && recoveryFlag == 1
        recoveryRT= [recoveryRT;IntanBehaviour.cueHitTrace(i).reactionTime];
    end
end

[pbc,~] = ranksum(coolRT,baselineRT)
[pcr,~] = ranksum(coolRT,recoveryRT)
[pbr,~] = ranksum(baselineRT,recoveryRT)
h=figure();plotBox3(baselineRT,coolRT,recoveryRT);
legend(strcat('pbc = ',num2str(pbc)),strcat('pcr = ',num2str(pcr)),strcat('pbr = ',num2str(pbr)))
xtix = {'Baseline','Cooled','Recovery'}; xtixloc = [1 2 3]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);set(gca,'TickDir','out','fontsize',14');
set(gca,'TickDir','out','fontsize',14');box off;
ylabel('Reaction Time (s)');
saveas(h,[savepath '\Figures\RTStatsCooling.png']);
saveas(h,[savepath '\Figures\RTStatsCooling.fig']);
%% Plotting wave speed 
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


[pbc,~] = ranksum(baselineSpeed',cooledSpeed')
[pcr,~] = ranksum(cooledSpeed',recoverySpeed')
[pbr,~] = ranksum(baselineSpeed',recoverySpeed')

temp = [baselineSpeed';cooledSpeed';recoverySpeed'];
templabels = cellstr([repmat('Baseline',size(baselineSpeed,2),1);repmat('Cooled  ',size(cooledSpeed,2),1);repmat('Recovery',size(recoverySpeed,2),1)]);
colors = [166/255 14/255 90/255;53/255 189/255 206/255;0.8500 0.3250 0.0980];
h=figure();violinplot(temp,templabels,'ShowData',true,'ShowWhiskers',false,'ShowBox',false,'MarkerSize',5,'ViolinColor',colors);
ylabel('Wave Speed (cm/s)'); title('Cooling M2');
% xtix = {'Baseline','Cooled','Recovery'}; xtixloc = [1 2 3]; set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);set(gca,'TickDir','out','fontsize',14');
set(gca,'TickDir','out','fontsize',14');box off;
legend(strcat('pbc = ',num2str(pbc)),strcat('pcr = ',num2str(pcr)),strcat('pbr = ',num2str(pbr)))
saveas(h,[savepath '\Figures\WaveSpeedStatsCooling.png']);
saveas(h,[savepath '\Figures\WaveSpeedStatsCooling.fig']);

% Waves speed vs time 
for i = 1:30
    st = (i-1)*100 + 1;
    sp = (i)*100 + 1;
    dirCombBaseline = horzcat(WavesBaseline.wavesHit(1:end).speed);
    evalPointsB = horzcat(WavesBaseline.wavesHit(1:end).evaluationPoints);
    WaveComb(i).Baseline = dirCombBaseline(evalPointsB >=st & evalPointsB <= sp);

    dirCombCool = horzcat(WavesCooling.wavesHit(1:end).speed);
    evalPointsC = horzcat(WavesCooling.wavesHit(1:end).evaluationPoints);
    WaveComb(i).Cool = dirCombCool(evalPointsC >=st & evalPointsC <= sp);

    dirCombRecovery = horzcat(WavesRecovery.wavesHit(1:end).speed);
    evalPointsR = horzcat(WavesRecovery.wavesHit(1:end).evaluationPoints);
    WaveComb(i).Recovery = dirCombRecovery(evalPointsR >=st & evalPointsR <= sp);
end

h=figure(); hold on;
avgVarBaseline = arrayfun(@(s) mean(s.Baseline,'all'), WaveComb);
maxVarBaseline = max(avgVarBaseline);
avgVarCool = arrayfun(@(s) mean(s.Cool,'all'), WaveComb);
maxVarCool = max(avgVarCool);
avgVarRecovery = arrayfun(@(s) mean(s.Recovery,'all'), WaveComb);
maxVarRecovery = max(avgVarRecovery);
abc = arrayfun(@(s) ranksum(s.Baseline,s.Cool) , WaveComb);
pVal = (abc<0.05);
significanceXbc = find(pVal)*100;
significanceYbc = 1.1*max([maxVarBaseline,maxVarCool,maxVarRecovery])*ones(1,length(significanceXbc));

acr = arrayfun(@(s) ranksum(s.Cool,s.Recovery) , WaveComb);
pVal = (acr<0.05);
significanceXcr = find(pVal)*100;
significanceYcr = 1.2*max([maxVarBaseline,maxVarCool,maxVarRecovery])*ones(1,length(significanceXcr));

abr = arrayfun(@(s) ranksum(s.Baseline,s.Recovery) , WaveComb);
pVal = (abr<0.05);
significanceXbr = find(pVal)*100;
significanceYbr = 1.3*max([maxVarBaseline,maxVarCool,maxVarRecovery])*ones(1,length(significanceXbr));

h1 = plot(100:100:3000,(avgVarBaseline),'Color', [166/255 14/255 90/255],'LineWidth',1.5);
h2 = plot(100:100:3000,(avgVarCool),'Color', [53/255 189/255 206/255],'LineWidth',1.5);
h3 = plot(100:100:3000,(avgVarRecovery),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
xline(1501,'--r','Cue');xlabel('Time (ms)'); ylabel('Wave speed (cm/s)');
plot(significanceXbc,significanceYbc,'r*');
plot(significanceXcr,significanceYcr,'g*');
plot(significanceXbr,significanceYbr,'b*');
legend([h1 h2 h3],'Baseline','Cool','Recovery','Location','best');
set(gca,'TickDir','out','fontsize',14')
title('M2 Cool');
saveas(h,[savepath '\Figures\WaveSpeedTimeCooling.png']);
saveas(h,[savepath '\Figures\WaveSpeedTimeCooling.fig']);

%% Get avg PGD for baseline, cooling and recovery 
PGD.Baseline.cueHits = vertcat(WavesBaseline.wavesHit.PGD);
PGD.Baseline.cueMiss = vertcat(WavesBaseline.wavesMiss.PGD);
PGD.Baseline.MIHits = vertcat(WavesBaseline.wavesMIHit.PGD);
PGD.Baseline.MIFA = vertcat(WavesBaseline.wavesMIFA.PGD);

PGD.Cooling.cueHits = vertcat(WavesCooling.wavesHit.PGD);
PGD.Cooling.cueMiss = vertcat(WavesCooling.wavesMiss.PGD);
PGD.Cooling.MIHits = vertcat(WavesCooling.wavesMIHit.PGD);
PGD.Cooling.MIFA = vertcat(WavesCooling.wavesMIFA.PGD);

PGD.Recovery.cueHits = vertcat(WavesRecovery.wavesHit.PGD);
PGD.Recovery.cueMiss = vertcat(WavesRecovery.wavesMiss.PGD);
PGD.Recovery.MIHits = vertcat(WavesRecovery.wavesMIHit.PGD);
PGD.Recovery.MIFA = vertcat(WavesRecovery.wavesMIFA.PGD);

PGD.PGDHits = vertcat(Waves.wavesHit.PGD);
PGD.PGDMiss = vertcat(Waves.wavesMiss.PGD); 
PGD.PGDMIHit = vertcat(Waves.wavesMIHit.PGD);
PGD.PGDMIFA = vertcat(Waves.wavesMIFA.PGD); 

% PGD vs cooling for hits 

h=figure();hold on;
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(mean(PGD.Baseline.cueHits,1),100,'sgolay',4),'Color', [166/255 14/255 90/255],'LineWidth',1.5);
plot(IntanBehaviourCooling.cueHitTrace(1).time,smooth(mean(PGD.Cooling.cueHits,1),100,'sgolay',4),'Color', [53/255 189/255 206/255],'LineWidth',1.5);
plot(IntanBehaviourRecovery.cueHitTrace(1).time,smooth(mean(PGD.Recovery.cueHits,1),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(mean(PGD.Baseline.cueHits,1)-(std(PGD.Baseline.cueHits,0,1)/sqrt(size(PGD.Baseline.cueHits,1))),100,'sgolay',4),'Color', [166/255 14/255 90/255 0.4],'LineWidth',0.5);
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(mean(PGD.Baseline.cueHits,1)+(std(PGD.Baseline.cueHits,0,1)/sqrt(size(PGD.Baseline.cueHits,1))),100,'sgolay',4),'Color', [166/255 14/255 90/255 0.4],'LineWidth',0.5);
plot(IntanBehaviourCooling.cueHitTrace(1).time,smooth(mean(PGD.Cooling.cueHits,1)-(std(PGD.Cooling.cueHits,0,1)/sqrt(size(PGD.Cooling.cueHits,1))),100,'sgolay',4),'Color', [53/255 189/255 206/255 0.4],'LineWidth',0.5);
plot(IntanBehaviourCooling.cueHitTrace(1).time,smooth(mean(PGD.Cooling.cueHits,1)+(std(PGD.Cooling.cueHits,0,1)/sqrt(size(PGD.Cooling.cueHits,1))),100,'sgolay',4),'Color', [53/255 189/255 206/255 0.4],'LineWidth',0.5);
plot(IntanBehaviourRecovery.cueHitTrace(1).time,smooth(mean(PGD.Recovery.cueHits,1)-(std(PGD.Recovery.cueHits,0,1)/sqrt(size(PGD.Recovery.cueHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviourRecovery.cueHitTrace(1).time,smooth(mean(PGD.Recovery.cueHits,1)+(std(PGD.Recovery.cueHits,0,1)/sqrt(size(PGD.Recovery.cueHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
% xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD) for Hits');box off;  legend('Baseline','Cooling','Recovery');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');
saveas(h,[savepath '\Figures\PGDHitsvsCooling.png']);
saveas(h,[savepath '\Figures\PGDHitsvsCooling.fig']);

h=figure();hold on;
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDHits,1),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDMiss,1),100,'sgolay',4),'Color', [0.7 0.7 0.7],'LineWidth',1.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDHits,1)-(std(PGD.PGDHits,0,1)/sqrt(size(PGD.PGDHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDHits,1)+(std(PGD.PGDHits,0,1)/sqrt(size(PGD.PGDHits,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDMiss,1)-(std(PGD.PGDMiss,0,1)/sqrt(size(PGD.PGDMiss,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDMiss,1)+(std(PGD.PGDMiss,0,1)/sqrt(size(PGD.PGDMiss,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Misses');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');
saveas(h,[savepath '\Figures\PGDHitsvsMiss.png']);
saveas(h,[savepath '\Figures\PGDHitsvsMiss.fig']);

h=figure();hold on;
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDMIHit,1),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980],'LineWidth',1.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDMIFA,1),100,'sgolay',4),'Color', [0.7 0.7 0.7],'LineWidth',1.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDMIHit,1)-(std(PGD.PGDMIHit,0,1)/sqrt(size(PGD.PGDMIHit,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDMIHit,1)+(std(PGD.PGDMIHit,0,1)/sqrt(size(PGD.PGDMIHit,1))),100,'sgolay',4),'Color', [0.8500 0.3250 0.0980 0.4],'LineWidth',0.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDMIFA,1)-(std(PGD.PGDMIFA,0,1)/sqrt(size(PGD.PGDMIFA,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
plot(IntanBehaviour.cueHitTrace(1).time,smooth(mean(PGD.PGDMIFA,1)+(std(PGD.PGDMIFA,0,1)/sqrt(size(PGD.PGDMIFA,1))),100,'sgolay',4),'Color', [0.7 0.7 0.7 0.4],'LineWidth',0.5);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','MI','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('MI Hits','MI FAs');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');
saveas(h,[savepath '\Figures\PGDHitsvsFA.png']);
saveas(h,[savepath '\Figures\PGDHitsvsFA.fig']);

sessionName = [savepath,'/','PGD.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"PGD","fpath","parameters","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",

%% Get avg PA for baseline, cooling and recovery 
PA.PABaseline = getPA(IntanBehaviourBaseline,0,1,0,parameters,0);
PA.PACooling = getPA(IntanBehaviourCooling,0,1,0,parameters,0);
PA.PARecovery = getPA(IntanBehaviourRecovery,0,1,0,parameters,0);

sessionName = [savepath,'/','PA.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"PA","fpath","parameters","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",

h=figure();
subplot(2,2,[1,2])
plot(IntanBehaviourBaseline.cueHitTrace(1).time,smooth(squeeze(mean(PA.PABaseline.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[166/255 14/255 90/255],'LineWidth',1.5); hold on;
plot(IntanBehaviourCooling.cueHitTrace(1).time,smooth(squeeze(mean(PA.PACooling.Hit,[1 2],'omitnan')),50,'sgolay',20),'Color',[53/255 189/255 206/255],'LineWidth',1.5);
ylabel("Phase Alignment"); xlabel("Time (s)");
xline(0,'--k','Cue','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xline((PA.PABaseline.PAPeakHit/parameters.Fs-parameters.windowBeforeCue),'--b','Peak','LabelVerticalAlignment','top','LabelHorizontalAlignment','left');
xlim([-0.5 1.5]);box off;legend('Baseline','Cool');legend('boxoff');set(gca,'TickDir','out','fontsize',14');
title("Phase Alignment : M2 -> M1 Cool");set(gca,'TickDir','out','fontsize',14');

subplot(2,2,3)
PABaselinePeak = PA.PABaseline.PAPeakHit;
xgp = arrayfun(@(s) s.xgp, IntanBehaviourBaseline.cueHitTrace, 'UniformOutput', false);
PABaselineAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PABaselinePeak)),parameters.rows*parameters.cols,1), xgp,'UniformOutput',false)),[],1));
histogram(PABaselineAngles,18,'FaceAlpha',0.7,'FaceColor',[166/255 14/255 90/255],'Normalization','probability','EdgeColor','none');hold on;
xgp = arrayfun(@(s) s.xgp, IntanBehaviourCooling.cueHitTrace, 'UniformOutput', false);
PACoolingAngles = rmmissing(reshape(cell2mat(cellfun(@(s) reshape(angle(s(:,:,PABaselinePeak)),parameters.rows*parameters.cols,1), xgp,'UniformOutput',false)),[],1));
histogram(PACoolingAngles,18,'FaceAlpha',0.75,'FaceColor',[53/255 189/255 206/255],'Normalization','probability','EdgeColor','none');box off;
xlabel('Peak PA Angle');ylabel('Probablitity');title('PA Angle at Peak');set(gca,'TickDir','out','fontsize',14');
subplot(2,2,4)
polarhistogram(PABaselineAngles,18,'FaceAlpha',0.7,'FaceColor',[166/255 14/255 90/255],'Normalization','probability','EdgeColor','none');hold on;
polarhistogram(PACoolingAngles,18,'FaceAlpha',0.75,'FaceColor',[53/255 189/255 206/255],'Normalization','probability','EdgeColor','none');box off;
title('PA Angle at Peak');set(gca,'TickDir','out','fontsize',14');

[p,~,~] = circ_kuipertest(PABaselineAngles, PACoolingAngles,60,0);
disp(['Peak PA angle p-val = ' num2str(p)]);

saveas(h,[savepath '\Figures\PAHitsvsCooling.png']);
saveas(h,[savepath '\Figures\PAHitsvsCooling.fig']);