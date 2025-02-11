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
addpath(genpath('Channel Maps'));
addpath(genpath('Dependancies'));
addpath(genpath('Pipelines'));
addpath(genpath('BatchProcessing'));
rmpath(genpath('Dependancies/MVGC1'));
%%  PreProcessing
load Grids30Ch_chanmap.mat;  % load the channel map for the IntanConcatenate function
% load GridsLowDenNeedle_chanmap.mat;  % load the channel map for the IntanConcatenate function
% load 64ChGrid_chanmap.mat;  % load the channel map for the IntanConcatenate function
% load 64ChM1M2Grid_chanmap.mat % load the channel map for 64ch Dual grid 
parameters.rows = 5;  % Number of rows of electrodes on the Grid
parameters.cols = 6;  % Number of colums of electrodes on the Grid
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;
parameters.windowBeforePull = 1.5; % in seconds
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds 
parameters.windowAfterCue = 1.5; % in seconds 
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated 
parameters.cool = 0; % 1 - cooling , 0 - no cooling 
parameters.opto = 1; % 1 - opto ON , 0 - opto OFF
parameters.xspacing = 0.12; % Grid spacing in mm between columns 
parameters.yspacing = 0.12; % Grid spacing in mm between rows
parameters.shank = 0; % 1 - if UCLA 64Ch single shank data is present
parameters.eOPN = 1; % 1 - if the opto trials were eOPN, 0 - otherwise

if parameters.shank == 1
    load UCLASingle64Ch_chanmap.mat; % load the channel map for the the shank data
    UCLAProbeMap = UCLAProbeMap + size(electrode_map,1);
    finalElectrodeMap = [electrode_map;UCLAProbeMap];
%     finalElectrodeMap = [UCLAProbeMap;electrode_map];
    parameters.nShank = size(UCLAProbeMap,1);
else
    finalElectrodeMap = electrode_map;
end

%% Getting path for arduino behaviour file
[enfile,enpath] = uigetfile('*.csv','Select csv file for Arduino behaviour file');
if isequal(enfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(enpath,enfile)]);
end

%% Getting path for impedance file
[impfile,imppath] = uigetfile('*.csv','Select csv file for Grid impedance');
if isequal(impfile,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(imppath,impfile)]);
end

%% Reading intan files 
IntanConcatenate
fpath = Intan.path; % where on disk do you want the analysis? ideally and SSD...

% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = 0:Ts:Intan.Tmax-Ts;

%% Removing bad channels from impedance values
% [Z,IntanBehaviour.goodChMap,IntanBehaviour.badChMap] = readImp(electrode_map,5e6);
[Z,IntanBehaviour.goodChMap,IntanBehaviour.badChMap] = readImp(imppath,impfile,finalElectrodeMap,5e6);
figure('Name','Impedance Test at 1kHz');boxchart(Z); xlabel('n = ' + string(size(Z,1)));ylabel('Impedance (in \Omega)');set(gca,'xticklabel',{[]})
% Intan.badChMap =[21,22];[1,2];[6,31];[5,10,21]; ;2;7];
%Intan = removeBadCh(Intan,Intan.badCh);
% IntanBehaviour.badChMap =[21,22];
% IntanBehaviour.badChMap =[31,32];
% IntanBehaviour.badChMap =[24,25,26];
% IntanBehaviour.badChMap =[29];
% IntanBehaviour.badChMap =[];
%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(Intan.allIntan,5000);
% LFP = bestLFP(LFP);
% LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
% LFPplot(LFP);
if parameters.shank == 1
    LFP.LFPprobe = reshape(LFP.LFP(size(electrode_map,1)+1:end,:),64,1,[]); % LFP from probe
    LFP.LFPprobe = LFP.LFPprobe(linearProbe,:,:);
end
% Creating LFPDataCube
LFP = createDataCube(LFP,parameters.rows,parameters.cols,IntanBehaviour.badChMap); % Creating datacube

%% Loading Lever Data 
plotOption = 0;
[Behaviour] = readLever(enpath,enfile,parameters,LFP.times,plotOption);

%% Reading behaviour data from Intan traces 
plotOption = 1;
IntanBehaviour = readLeverIntan(parameters,LFP.times,Intan.analog_adc_data,Intan.dig_in_data,Behaviour,plotOption);

%% Generalized Phase 
[parameters.X,parameters.Y] = meshgrid( 1:parameters.cols, 1:parameters.rows );
IntanBehaviour = addLFPToBehaviour(IntanBehaviour,LFP,parameters);

%% Wave detection in trials
nShuffle = 50;
threshold = 99.73; % zscore of 3
fraction = 0.1;
parameters.rhoThres = getRhoThreshold(IntanBehaviour.cueHitTrace,IntanBehaviour.cueMissTrace,parameters,nShuffle,threshold,fraction);
% parameters.rhoThres = 0.75;

disp('Wave Detection for wide band ...')
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
% Analysis code for wave detection on beta,theta,gamma band
% multibandWaveAnalysis;

%%
%% Saving paramters, path, IntanBehaviour to bin file 
savepath = uigetdir(path);
sessionName = [savepath,'/','61590eOPNDay2Waves.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"IntanBehaviour","fpath","parameters","Waves","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",

if parameters.opto == 1
    [IntanBehaviourBaseline,IntanBehaviourOpto, WavesBaseline, WavesOpto] = separateOptoTrials(IntanBehaviour,parameters,Waves);
    clear('IntanBehaviour','Waves');
    save(sessionName,"IntanBehaviourBaseline","IntanBehaviourOpto","fpath","parameters","WavesBaseline","WavesOpto","-v7.3");
end

%% Wave detecion for entire time 
parameters.rhoThres= 0.75;
allwaves.LFPIndex = (1:1:size(LFP.LFP,2))';
xf{1,1} = LFP.xf;
xgp{1,1} = LFP.xgp;
wt{1,1} = LFP.wt;
Wavesall = detectWaves(xf,xgp,wt,allwaves,parameters,0.6);
WaveStatsSingle(Wavesall,parameters,1,size(LFP.LFP,2)/20);

%% PLotting to check visually
% trialPlot = 13;
% plot_wave_examples( LFP.xf(:,:,IntanBehaviour.cueHitTrace(trialPlot).LFPIndex(1):IntanBehaviour.cueHitTrace(trialPlot).LFPIndex(end)), ...
%     options, trialPlot, Waves.wavesHit,rhoThres);

% trialPlot = 13;
% plot_wave_examples( LFP.xf(:,:,IntanBehaviour.missTrace(trialPlot).LFPIndex(1):IntanBehaviour.missTrace(trialPlot).LFPIndex(end)), ...
%     options, trialPlot, Waves.wavesFA,rhoThres);

%% Waves accross trials 
plotOption = 1;
[WaveStats(1)] = getWaveStats(Waves.wavesHit,Waves.wavesMiss,parameters,plotOption);
[WaveStats(2)] = getWaveStats(thetaWaves.wavesHit,thetaWaves.wavesMiss,parameters,plotOption);
[WaveStats(3)] = getWaveStats(betaWaves.wavesHit,betaWaves.wavesMiss,parameters,plotOption);
[WaveStats(4)] = getWaveStats(gammaWaves.wavesHit,gammaWaves.wavesMiss,parameters,plotOption);

plotOption = 1;
[WaveStatsReward(1)] = getWaveStats(Waves.wavesHitReward,Waves.wavesFA,parameters,plotOption);

plotOption = 1;
[WaveStatsFA(1)] = getWaveStats(Waves.wavesMIHit,Waves.wavesMIFA,parameters,plotOption);

[WaveStats2(1)] = getInitRewardStats(Waves,parameters,plotOption);
[WaveStats2(2)] = getInitRewardStats(thetaWaves,parameters,plotOption);
[WaveStats2(3)] = getInitRewardStats(betaWaves,parameters,plotOption);
[WaveStats2(4)] = getInitRewardStats(gammaWaves,parameters,plotOption);

[WaveStats(1)] = getWaveStats(betaWaves.wavesHit,betaWaves.wavesMiss,parameters,plotOption);
[WaveStats(1)] = getWaveStats(gammaWaves.wavesHit,gammaWaves.wavesMiss,parameters,plotOption);

%% Combining  multiple Intanbehaviour structs from multiple sessions
% combIntanBehaviour = horzcat(IntanBehaviourBaseline, IntanBehaviourCool);
% IntanBehaviour.cueHitTrace = horzcat(combIntanBehaviour(1:end).cueHitTrace);
% IntanBehaviour.cueMissTrace = horzcat(combIntanBehaviour(1:end).cueMissTrace);
% IntanBehaviour.hitTrace = horzcat(combIntanBehaviour(1:end).hitTrace);
% IntanBehaviour.missTrace = horzcat(combIntanBehaviour(1:end).missTrace);
% IntanBehaviour.MIHitTrace = horzcat(combIntanBehaviour(1:end).MIHitTrace);
% IntanBehaviour.MIFATrace = horzcat(combIntanBehaviour(1:end).MIFATrace);
% IntanBehaviour.reactionTime = horzcat(combIntanBehaviour(1:end).reactionTime);
% 
% clear combIntanBehaviour %IntanBehaviour1 IntanBehaviour2;
%% Power Spectrum during task across channels 
PSDAnalyis;

%% Wavelet spectrogram
waveletSpectrogramAnalysis;

%% Mutual Information
z_score = 0;
nIterrate = 10;
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
plot(IntanBehaviour.cueHitTrace(1).time,smooth(PGD.avgPGDHit,50,'sgolay',20),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.cueMissTrace(1).time,smooth(PGD.avgPGDMiss,50,'sgolay',20),'-k','LineWidth',1);
% plot(IntanBehaviour.cueHitTrace(1).time,PGD.avgPGDFA,'-b','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','Misses');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');


figure(); hold on;
plot(IntanBehaviour.hitTrace(1).time,PGD.avgPGDHitReward,'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.missTrace(1).time,PGD.avgPGDFA,'-k','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','Threshold','LabelVerticalAlignment','top');
xline(-mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Cue Time','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD)');box off;  legend('Hits','FAs');%,'False Alarms');

figure(); hold on;
plot(IntanBehaviour.MIHitTrace(1).time,smooth(PGD.avgPGDMIHit,50,'sgolay',20),'-r','LineWidth',1.2); hold on;
plot(IntanBehaviour.MIFATrace(1).time,smooth(PGD.avgPGDMIFA,50,'sgolay',20),'-k','LineWidth',1);
ylabel("PGD"); xlabel("Time (s)");
xline(0,'--r','MI','LabelVerticalAlignment','top');
xline(-mean(IntanBehaviour.reactionTime,'all'),'--m','RT','LabelVerticalAlignment','top');
title('Trial Averaged Phase Gradient  Directionality (PGD) - Motion Ininitiation');box off;  legend('Hits','FAs');%,'False Alarms');
xlim([-0.5 1.5]);set(gca,'TickDir','out','fontsize',14');

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
z_score = 1;
nIterrate = 100;
PAGrid = getPA(IntanBehaviour,z_score,nIterrate,1,parameters,0);
PAProbe = getPA(IntanBehaviour,z_score,nIterrate,1,parameters,1);

% PA = PAGrid;
% PA= PAProbe;
%% Percent Phase Locking
z_score = 1;
nIterrate = 10;
PPLGrid = getPPL(IntanBehaviour,z_score,nIterrate,1,parameters,0);
PPLShank = getPPL(IntanBehaviour,z_score,nIterrate,1,parameters,1);

% PPL = PPLGrid;
% PPL = PPLShank;
%% Average LFP for Hits and Misses 
LFPHit = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.cueHitTrace,"UniformOutput",false);
LFPMiss = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.cueMissTrace,"UniformOutput",false);
LFPHitProbe = arrayfun(@(s) reshape(s.xfProbe,[],size(s.xfProbe,3)), IntanBehaviour.cueHitTrace,"UniformOutput",false);
LFPMissProbe = arrayfun(@(s) reshape(s.xfProbe,[],size(s.xfProbe,3)), IntanBehaviour.cueMissTrace,"UniformOutput",false);
LFPHitReward = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.hitTrace,"UniformOutput",false);
LFPFA = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.missTrace,"UniformOutput",false);
LFPMIHit = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.MIHitTrace,"UniformOutput",false);
LFPMIFA = arrayfun(@(s) reshape(s.xf,[],size(s.xf,3)), IntanBehaviour.MIFATrace,"UniformOutput",false);

avgLFPHit = zeros(parameters.rows*parameters.cols,size(LFPHit{1,1},2));
avgLFPMiss = zeros(parameters.rows*parameters.cols,size(LFPMiss{1,1},2));
avgLFPHitProbe = zeros(size(LFPHit{1,1},1),size(LFPHit{1,1},2));
avgLFPMissProbe = zeros(size(LFPMiss{1,1},1),size(LFPMiss{1,1},2));
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

a = cell2struct(LFPHitProbe,'lfp',1);
for i=1:(parameters.nShank)
    avgLFPHitProbe(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
end

% Getting trial averaged LFP for each channel for Miss
a = cell2struct(LFPMissProbe,'lfp',1);
for i=1:(parameters.nShank)
    avgLFPMissProbe(i,:) = mean(cell2mat(arrayfun(@(s) s.lfp(i,:),a, 'UniformOutput',false)),1,'omitnan');
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

avgLFPHit = removeNaNRows(peakSort2DArray(avgLFPHit,'descend',2));
avgLFPMiss = removeNaNRows(peakSort2DArray(avgLFPMiss,'descend',2));

figure(); % Top half is hits and bottom half is misses
title("Trial Average LFP for Hits and Misses")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:size([avgLFPHit;avgLFPMiss],1),[avgLFPHit;avgLFPMiss]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','Cue','LabelVerticalAlignment','top');
% yline(32.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);

figure(); % Top half is hits and bottom half is misses for probes
title("Trial Average LFP for Hits and Misses")
imagesc(IntanBehaviour.cueHitTrace(1).time,1:2*parameters.nShank,[avgLFPHitProbe;avgLFPMissProbe]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','Cue','LabelVerticalAlignment','top');
yline(parameters.nShank+0.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);

avgLFPHitReward = removeNaNRows(peakSort2DArray(avgLFPHitReward,'descend',2));
avgLFPFA = removeNaNRows(peakSort2DArray(avgLFPFA,'descend',2));

figure(); % Top half is hits and bottom half is FA
title("Trial Average LFP for Hits and FA")
imagesc(IntanBehaviour.hitTrace(1).time,1:size([avgLFPHitReward;avgLFPFA],1),[avgLFPHitReward;avgLFPFA]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','Reward','LabelVerticalAlignment','top');
yline(32.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);


avgLFPMIHit = removeNaNRows(peakSort2DArray(avgLFPMIHit,'descend',2));
avgLFPMIFA = removeNaNRows(peakSort2DArray(avgLFPMIFA,'descend',2));

figure(); % Top half is hits MI and bottom half is FA MI
title("Trial Average LFP for Hits and FA- Motion Initiation")
imagesc(IntanBehaviour.MIHitTrace(1).time,1:size([avgLFPMIHit;avgLFPMIFA],1),[avgLFPMIHit;avgLFPMIFA]); colormap(jet);
ylabel("Electrodes");xlabel("Time (s)"); 
h = colorbar; h.Label.String = 'Amplitude (uV)';
xline(0,'-k','MI','LabelVerticalAlignment','top');
yline(32.5,'-k');caxis([-30 30]);
% yline(64.5,'-k');caxis([-20 20]);
%% Average PGD accross frequency bands
[PGDfreqHit,PGDfreqMiss] = getPGDFreqBand(LFP.LFPdatacube,IntanBehaviour,0,parameters);

figure();
PGDFreq = [5:5:100];
plot(PGDFreq,PGDfreqHit); hold on;
plot(PGDFreq,PGDfreqMiss);
xlabel('Frequency'); ylabel('Phase Gradient Directionality'); ylim([0.35 0.7]); box off;


%% Average PA accross frequency bands 

nMax = 45;
PAfreqHit = zeros(nMax+1,parameters.rows*parameters.cols,(parameters.windowAfterCue+parameters.windowBeforeCue)*parameters.Fs+1);
PAfreqMiss = zeros(nMax+1,parameters.rows*parameters.cols,(parameters.windowAfterCue+parameters.windowBeforeCue)*parameters.Fs+1);
xgpcellHit = cell(1,size(IntanBehaviour.cueHitTrace,2));
xgpcellMiss = cell(1,size(IntanBehaviour.cueMissTrace,2));

for i=1:nMax
    xf = bandpass_filter(LFP.LFPdatacube,i+4,(i+5),4,1000);
    [xgp, ~] = generalized_phase(xf,1000,0);
    for jj=1:size(IntanBehaviour.cueHitTrace,2)
       xgpcellHit{1,jj} = reshape(xgp(:,:,IntanBehaviour.cueHitTrace(jj).LFPIndex(1):IntanBehaviour.cueHitTrace(jj).LFPIndex(end)),parameters.cols*parameters.rows,[]);
    end
    for jj=1:size(IntanBehaviour.cueMissTrace,2)
       xgpcellMiss{1,jj} = reshape(xgp(:,:,IntanBehaviour.cueMissTrace(jj).LFPIndex(1):IntanBehaviour.cueMissTrace(jj).LFPIndex(end)),parameters.cols*parameters.rows,[]);
    end
    PAfreqHit(i,:,:) = calPhaseAlignment(xgpcellHit,parameters);
    PAfreqMiss(i,:,:) = calPhaseAlignment(xgpcellMiss,parameters);
    i
end

freqPA = 5:1:nMax+5;
figure();
subplot(1,2,1)
imagesc(IntanBehaviour.cueHitTrace(1).time,freqPA,squeeze(mean(PAfreqHit,2,'omitnan')));
axis xy;
subplot(1,2,2)
imagesc(IntanBehaviour.cueMissTrace(1).time,freqPA,squeeze(mean(PAfreqMiss,2,'omitnan')));
axis xy;


freqPA = 5:1:nMax+5;
figure();
plot(IntanBehaviour.cueHitTrace(1).time,squeeze(mean(PAfreqHit(10:16,:,:),[1 2],'omitnan')));
plot(IntanBehaviour.cueMissTrace(1).time,squeeze(mean(PAfreqMiss(10:16,:,:),[1 2],'omitnan')));

%% Beta Burst detection using Hilbert Amplitude 

[GammaEventHit] = detectBetaEvent(LFP.xgpgamma,IntanBehaviour.MIHitTrace,parameters);
[GammaEventMiss] = detectBetaEvent(LFP.xgpgamma,IntanBehaviour.MIFATrace,parameters);

[BetaEventHit] = detectBetaEvent(LFP.xgpbeta,IntanBehaviour.MIHitTrace,parameters);
[BetaEventMiss] = detectBetaEvent(LFP.xgpbeta,IntanBehaviour.MIFATrace,parameters);

%% Beta event detection 
avgBetaband = mean(LFP.xfbeta,1);
window = Encoder.trialTime(:,3:4);
windowStop = Encoder.trialTimeStop(:,3:4);

betatrials = zeros(Encoder.nTrig,LFP.Fs*(windowAfterTrig+windowBeforeTrig));
betatrialsStop = zeros(Encoder.nTrigStop,LFP.Fs*(windowAfterTrig+windowBeforeTrig));
for i=1:Encoder.nTrig
   betatrials(i,:) = avgBetaband(window(i,1):window(i,2));
end
for i=1:Encoder.nTrigStop
   betatrialsStop(i,:) = avgBetaband(windowStop(i,1):windowStop(i,2));
end
avgbetaGroup = groupBetaBurstDetection(LFP,betatrials',window,LFP.Fs);
avgbetaGroupStop = groupBetaBurstDetection(LFP,betatrialsStop',windowStop,LFP.Fs);

% Calculating for each electrode
betatrials = zeros(Encoder.nTrig,LFP.Fs*(windowAfterTrig+windowBeforeTrig));
betatrialsStop = zeros(Encoder.nTrigStop,LFP.Fs*(windowAfterTrig+windowBeforeTrig));
for i=1:parameters.rows*parameters.cols
    for j=1:Encoder.nTrig
        betatrials(j,:) = LFP.beta_band(i,window(j,1):window(j,2));
    end
    for j=1:Encoder.nTrigStop
        betatrialsStop(j,:) = LFP.beta_band(i,windowStop(j,1):windowStop(j,2));
    end
    betaEvents(i).betagroup = groupBetaBurstDetection(LFP,betatrials',window,LFP.Fs);
    betaEventsStop(i).betagroup = groupBetaBurstDetection(LFP,betatrialsStop',windowStop,LFP.Fs);
end

% Number of beta events per trial
nbetaevents = zeros(Encoder.nTrig,1);
nbetaeventsStop = zeros(Encoder.nTrigStop,1);
for i=1:parameters.rows*parameters.cols
    nbetaevents = nbetaevents + betaEvents(i).betagroup.betaBurst.NumDetectedBeta;
    nbetaeventsStop = nbetaeventsStop + betaEventsStop(i).betagroup.betaBurst.NumDetectedBeta;
end
figure,plot(nbetaevents);hold on;plot(nbetaeventsStop);

% Number of beta events on each electrode 
nbetaperElec = zeros(parameters.cols,parameters.rows);
for i=1:parameters.rows*parameters.cols
    a = chToGrid(i,parameters);
    nbetaperElec(a(1),a(2)) = mean(betaEvents(i).betagroup.betaBurst.NumDetectedBeta,"all");
end

figure(); imagesc(nbetaperElec');%set(gca,'YDir','normal');
title('Spatial map of beta events accross all trials'); colorbar;

%% PLotting LFP 
figure();
trialno = 1;
for i=1:32
    subplot(parameters.rows,parameters.cols,i);
    if ismember(i,Intan.badChMap), continue; end
    plot(LFP.times(Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4)),squeeze(goodTrial.xf(floor((i-1)/parameters.cols)+1,mod(i-1,parameters.cols)+1,:))');xline(0,'-r');
end
