% clear; clc; 
% close all;
format compact;
% set(0,'DefaultFigureWindowStyle','normal')

addpath(genpath('main'));
addpath(genpath('chronux'));
addpath(genpath('Kilosort'));
addpath(genpath('npy-matlab'));
addpath(genpath('spikes-master'));
addpath(genpath('PreProcessing'));
addpath(genpath('Plotting'));
addpath(genpath('Analysis'));
addpath(genpath('Dependancies'));
rmpath(genpath('Dependancies/MVGC1'));

%% Wave detection 
nShuffle = 1000;
threshold = 99;

% Wave detection for wide band
disp('Wave Detection for wide band ...')
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
parameters.rhoThres = getRhoThreshold(xgp,IntanBehaviour.hitTrace,parameters,nShuffle,trialno,threshold);
if isfield(IntanBehaviour,'cueHitTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
    Waves.wavesHit = detectWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,parameters.rhoThres);
end
% Waves.wavesHit = detectPlanarWaves(xf,xgp,wt,IntanBehaviour.cueHitTrace,parameters,0.5);
if isfield(IntanBehaviour,'cueMissTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
    Waves.wavesMiss = detectWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,parameters.rhoThres);
%     Waves.wavesMiss = detectPlanarWaves(xf,xgp,wt,IntanBehaviour.cueMissTrace,parameters,0.5);
end
if isfield(IntanBehaviour,'MIHitTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
    Waves.wavesMIHit = detectWaves(xf,xgp,wt,IntanBehaviour.MIHitTrace,parameters,parameters.rhoThres);
end
if isfield(IntanBehaviour,'MIFATrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.MIFATrace, 'UniformOutput', false);
    Waves.wavesMIFA = detectWaves(xf,xgp,wt,IntanBehaviour.MIFATrace,parameters,parameters.rhoThres);
end
if isfield(IntanBehaviour,'missTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.missTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.missTrace, 'UniformOutput', false);
    Waves.wavesFA = detectWaves(xf,xgp,wt,IntanBehaviour.missTrace,parameters,parameters.rhoThres);
end
if isfield(IntanBehaviour,'hitTrace')
    xf = arrayfun(@(s) s.xf, IntanBehaviour.hitTrace, 'UniformOutput', false);
    xgp = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
    wt = arrayfun(@(s) s.wt, IntanBehaviour.hitTrace, 'UniformOutput', false);
    Waves.wavesHitReward = detectWaves(xf,xgp,wt,IntanBehaviour.hitTrace,parameters,parameters.rhoThres);
end

%% PLotting to check visually
% trialPlot = 48;
% plot_wave_examples( LFP.xf(:,:,IntanBehaviour.cueHitTrace(trialPlot).LFPIndex(1):IntanBehaviour.cueHitTrace(trialPlot).LFPIndex(end)), ...
%     options, trialPlot, Waves.wavesHit,rhoThres);

% trialPlot = 12;
% plot_wave_examples( LFP.xf(:,:,IntanBehaviour.missTrace(trialPlot).LFPIndex(1):IntanBehaviour.missTrace(trialPlot).LFPIndex(end)), ...
%     options, trialPlot, Waves.wavesMiss,rhoThres);

%% Waves accross trials 
plotOption = 1;
[WaveStats(1)] = getWaveStats(Waves.wavesHit,Waves.wavesMiss,parameters,plotOption);

plotOption = 1;
[WaveStatsReward(1)] = getWaveStats(Waves.wavesHitReward,Waves.wavesFA,parameters,plotOption);

plotOption = 1;
[WaveStatsFA(1)] = getWaveStats(Waves.wavesMIHit,Waves.wavesMIFA,parameters,plotOption);


%% Mutual Information
z_score = 0;
nIterrate = 100;
MI = getMI(IntanBehaviour,z_score,nIterrate,1,parameters);

% For amplitude
xgpHit = arrayfun(@(s) abs(s.xgp), IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) abs(s.xgp), IntanBehaviour.cueMissTrace, 'UniformOutput', false);
[MI.Amp] = getMutualInformation(xgpHit,xgpMiss,parameters);

figure();
title("Mututal Information across all electrodes - Amplitude")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,peakSort2DArray(reshape(MI.Amp,[],size(MI.Amp,3)),'descend',2)); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Information (bits)';
xline(0,'-w','Cue','LabelVerticalAlignment','top');


%% Average PGD 
PGD.avgPGDHit = mean(vertcat(Waves.wavesHit.PGD),1);
PGD.avgPGDMiss = mean(vertcat(Waves.wavesMiss.PGD),1);
PGD.avgPGDHitReward = mean(vertcat(Waves.wavesHitReward.PGD),1);
PGD.avgPGDFA = mean(vertcat(Waves.wavesFA.PGD),1);
PGD.avgPGDMIHit = mean(vertcat(Waves.wavesMIHit.PGD),1);
PGD.avgPGDMIFA = mean(vertcat(Waves.wavesMIFA.PGD),1);

figure(); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDHit,'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDMiss,'-k','LineWidth',1);
% plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDFA,'-b','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Misses');%,'False Alarms');

figure(); hold on;
plot(IntanBehaviour.hitTrace(1).time,PGD.avgPGDHitReward,'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.missTrace(1).time,PGD.avgPGDFA,'-k','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Threshold','LabelVerticalAlignment','top');
xline(-mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','FAs');%,'False Alarms');

figure(); hold on;
plot(IntanBehaviour.MIHitTrace(1).time,PGD.avgPGDMIHit,'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.MIFATrace(1).time,PGD.avgPGDMIFA,'-k','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','MI','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD) - Motion Ininitiation');box off;  legend('Hits','FAs');%,'False Alarms');

xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
[PGD.avgPGDHit,PGD.avgPGDHitNull] = getAvgPGD(xgp,Waves.wavesHit,IntanBehaviour.cueHitTrace,1,parameters);

figure(); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDHit,'-k','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(PGD.avgPGDHitNull),'-r','LineWidth',1.2);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Shuffled');

%% Cross-Trial Phase Alignment
z_score = 0;
nIterrate = 2000;
PA = getPA(IntanBehaviour,z_score,nIterrate,1,parameters);

%% Percent Phase Locking
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
[PPLHit] = getPPL(xgp,parameters);
% xgp = arrayfun(@(s) shuffle3DMatrix(s.xgp,3), IntanBehaviour.cueHitTrace, 'UniformOutput', false);
% [PPLHitShuffle] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
[PPLMiss] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
[PPLHitReward] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
[PPLFA] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIHitTrace, 'UniformOutput', false);
[PPLMIHit] = getPPL(xgp,parameters);
xgp = arrayfun(@(s) s.xgp, IntanBehaviour.MIFATrace, 'UniformOutput', false);
[PPLMIFA] = getPPL(xgp,parameters);

figure();
subplot(2,1,1);
title("Percentage Phase across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PPLHit,[],size(PPLHit,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--w','Avg. Reaction Time','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Percentage Phase across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PPLMiss,[],size(PPLMiss,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');

figure();
subplot(2,1,1);
title("Percentage Phase across all electrodes - Hits")
imagesc(IntanBehaviour.hitTrace(1).time,1:32,reshape(PPLHitReward,[],size(PPLHitReward,3))); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Threshold','LabelVerticalAlignment','top');
yyaxis right; box off;
plot(IntanBehaviour.hitTrace(1).time,squeeze(mean(PPLHitReward,[1 2],'omitnan')),'-r','Linewidth',0.8);
% xline(-mean(IntanBehaviour.reactionTime,'all'),'--w','Avg. Cue Time','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Percentage Phase across all electrodes - FA")
imagesc(IntanBehaviour.missTrace(1).time,1:32,reshape(PPLFA,[],size(PPLFA,3))); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Threshold','LabelVerticalAlignment','top');
yyaxis right; box off;
plot(IntanBehaviour.missTrace(1).time,squeeze(mean(PPLFA,[1 2],'omitnan')),'-r','Linewidth',0.8);

figure();
subplot(2,1,1);
title("Percentage Phase across all electrodes - Hits MI")
imagesc(IntanBehaviour.MIHitTrace(1).time,1:32,reshape(PPLMIHit,[],size(PPLMIHit,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','MI','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Percentage Phase across all electrodes - Misses MI")
imagesc(IntanBehaviour.MIFATrace(1).time,1:32,reshape(PPLMIFA,[],size(PPLMIFA,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','MI','LabelVerticalAlignment','top');


figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLHit,[1 2])),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLMiss,[1 2])),'-k','LineWidth',0.2);
ylabel("Percentage Phase Locking"); xlabel("Time (s)");
xline(0,'--k','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--k','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Percentage Phase Locking for hits');box off;legend('Hits','Misses');%,'Misses');

figure();
plot(IntanBehaviour.MIHitTrace(1).time,squeeze(nanmean(PPLMIHit,[1 2])),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.MIFATrace(1).time,squeeze(nanmean(PPLMIFA,[1 2])),'-k','LineWidth',0.2);
ylabel("Percentage Phase Locking"); xlabel("Time (s)");
xline(0,'--k','MI','LabelVerticalAlignment','top');
title('Percentage Phase Locking for hits');box off;legend('Hits','FAs');%,'Misses');

% z-scoring 
nIterrate = 200;
xgpHit = arrayfun(@(s) s.xgp, IntanBehaviour.cueHitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) s.xgp, IntanBehaviour.cueMissTrace, 'UniformOutput', false);
xgpComb = [xgpHit xgpMiss];
nHit = size(xgpHit,2);
nMiss = size(xgpMiss,2);
nTot = nHit + nMiss;
nullDistHit = zeros(parameters.rows,parameters.cols,size(PPLHit,3),nIterrate);
nullDistMiss = zeros(parameters.rows,parameters.cols,size(PPLMiss,3),nIterrate);
for j=1:nIterrate
    randIndex = randperm(nTot);
    xgpHitRand = xgpComb(randIndex(1:nHit));
    xgpMissRand = xgpComb(randIndex(nHit+1:end));
    nullDistHit(:,:,:,j) = getPPL(xgpHitRand,parameters);
    nullDistMiss(:,:,:,j) = getPPL(xgpMissRand,parameters);
    j
end
muHit = mean(nullDistHit,4); % Mean of the null distribution
sigmaHit = std(nullDistHit,0,4); % Standard deviation of null distribution
muMiss = mean(nullDistMiss,4); % Mean of the null distribution
sigmaMiss = std(nullDistMiss,0,4); % Standard deviation of null distribution

PPLHitz = (PPLHit-muHit)./sigmaHit;
PPLMissz = (PPLMiss-muMiss)./sigmaMiss;

figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLHitz,[1 2])),'-r','LineWidth',1.2); hold on;
% plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PAMiss,[1 2])),'-k','LineWidth',1);
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(nanmean(PPLMissz,[1 2])),'-k','LineWidth',1);
ylabel("z-score"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('z-scored Phase Alignment for Hits');box off; legend('Hits','Miss');

figure();
subplot(2,1,1);
title("Percentage Phase across all electrodes - Hits")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:32,reshape(PPLHitz,[],size(PPLHitz,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');
subplot(2,1,2);
title("Percentage Phase across all electrodes - Misses")
imagesc(IntanBehaviour.cueMissTrace(1).time,1:32,reshape(PPLMissz,[],size(PPLMissz,3))); colormap(hot);
ylabel("Electrodes");xlabel("Time (s)");
xline(0,'-w','Cue','LabelVerticalAlignment','top');

nIterrate = 200;
xgpHit = arrayfun(@(s) s.xgp, IntanBehaviour.hitTrace, 'UniformOutput', false);
xgpMiss = arrayfun(@(s) s.xgp, IntanBehaviour.missTrace, 'UniformOutput', false);
xgpComb = [xgpHit xgpMiss];
nHit = size(xgpHit,2);
nMiss = size(xgpMiss,2);
nTot = nHit + nMiss;
nullDistHit = zeros(parameters.rows,parameters.cols,size(PPLHitReward,3),nIterrate);
nullDistMiss = zeros(parameters.rows,parameters.cols,size(PPLFA,3),nIterrate);
for j=1:nIterrate
    randIndex = randperm(nTot);
    xgpHitRand = xgpComb(randIndex(1:nHit));
    xgpMissRand = xgpComb(randIndex(nHit+1:end));
    nullDistHit(:,:,:,j) = getPPL(xgpHitRand,parameters);
    nullDistMiss(:,:,:,j) = getPPL(xgpMissRand,parameters);
    j
end
muHitReward = mean(nullDistHit,4); % Mean of the null distribution
sigmaHitReward = std(nullDistHit,0,4); % Standard deviation of null distribution
muFA = mean(nullDistMiss,4); % Mean of the null distribution
sigmaFA = std(nullDistMiss,0,4); % Standard deviation of null distribution

PPLHitRewardz = (PPLHit-muHitReward)./sigmaHitReward;
PPLFAz = (PPLMiss-muFA)./sigmaFA;


%% Average LFP for Hits and Misses 
LFPHit = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.cueHitTrace,"UniformOutput",false);
LFPMiss = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.cueMissTrace,"UniformOutput",false);
LFPHitReward = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.hitTrace,"UniformOutput",false);
LFPFA = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.missTrace,"UniformOutput",false);
LFPMIHit = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.MIHitTrace,"UniformOutput",false);
LFPMIFA = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.MIFATrace,"UniformOutput",false);
avgLFPHit = zeros(parameters.rows*parameters.cols,size(LFPHit{1,1},2));
avgLFPMiss = zeros(parameters.rows*parameters.cols,size(LFPMiss{1,1},2));
avgLFPHitReward = zeros(parameters.rows*parameters.cols,size(LFPHitReward{1,1},2));
avgLFPFA = zeros(parameters.rows*parameters.cols,size(LFPFA{1,1},2));
avgLFPMIHit = zeros(parameters.rows*parameters.cols,size(LFPMIHit{1,1},2));
avgLFPMIFA = zeros(parameters.rows*parameters.cols,size(LFPMIFA{1,1},2));

% Getting trial averaged LFP for each channel for Hits
a = cell2struct(LFPHit,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPHit(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for Miss
a = cell2struct(LFPMiss,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPMiss(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for Hits Rewards
a = cell2struct(LFPHitReward,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPHitReward(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for False alarms
a = cell2struct(LFPFA,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPFA(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for Hits MI
a = cell2struct(LFPMIHit,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPMIHit(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for False alarms MI
a = cell2struct(LFPMIFA,'lfp',1);
for i=1:(parameters.rows*parameters.cols)
    avgLFPMIFA(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

figure(); % Top half is hits and bottom half is misses
title("Trial Average LFP for Hits and Misses")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:64,[avgLFPHit;avgLFPMiss]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','Cue','LabelVerticalAlignment','top');
yline(32.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);

figure(); % Top half is hits and bottom half is FA
title("Trial Average LFP for Hits and FA")
imagesc(IntanBehaviour.hitTrace(1).time,1:64,[avgLFPHitReward;avgLFPFA]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','Reward','LabelVerticalAlignment','top');
yline(32.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);


figure(); % Top half is hits MI and bottom half is FA MI
title("Trial Average LFP for Hits and FA- Motion Initiation")
imagesc(IntanBehaviour.MIHitTrace(1).time,1:64,[avgLFPMIHit;avgLFPMIFA]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','MI','LabelVerticalAlignment','top');
yline(32.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);
