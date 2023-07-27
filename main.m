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

%% Loading previous workspace
% loadPreviousWS

%%  PreProcessing
load GridsLowDenNeedle_chanmap.mat;  % load the channel map for the IntanConcatenate function
parameters.rows = 8;  % Number of rows of electrodes on the Grid
parameters.cols = 4;  % Number of colums of electrodes on the Grid
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;
parameters.windowBeforePull = 1; % in seconds
parameters.windowAfterPull = 1; % in seconds
parameters.windowBeforeCue = 0.5; % in seconds 
parameters.windowAfterCue = 1.5; % in seconds 
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated 
parameters.opto = 1; % 1 - opto ON , 0 - opto OFF

IntanConcatenate
fpath = Intan.path; % where on disk do you want the analysis? ideally and SSD...

% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = 0:Ts:Intan.Tmax-Ts;

%% Removing bad channels from impedance values
[Z,Intan.goodChMap,Intan.badChMap] = readImp(electrode_map,100e6);
figure('Name','Impedance Test at 1kHz');boxchart(Z); xlabel('n = ' + string(size(Z,1)));ylabel('Impedance (in \Omega)');set(gca,'xticklabel',{[]})
%Intan = removeBadCh(Intan,Intan.badCh);

%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(Intan.allIntan,5000);
% LFP = bestLFP(LFP);
% LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
% LFPplot(LFP);
LFP = createDataCube(LFP,parameters.rows,parameters.cols,Intan.goodChMap); % Creating datacube

%% Loading Lever Data 
[Behaviour] = readLever(parameters,LFP.times);
figure('Name','Lever Trace');plot(Behaviour.time,Behaviour.leverTrace,'LineWidth',1.5);ylim([-5 50]);xlabel('Time (in s)');ylabel('Lever Position in mV');yline(23);
xline(squeeze(Behaviour.hit(:,2)),'-.b',cellstr(num2str((1:1:Behaviour.nHit)')),'LabelVerticalAlignment','top');
xline(squeeze(Behaviour.miss(:,2)),'-.r',cellstr(num2str((1:1:Behaviour.nMiss)')),'LabelVerticalAlignment','bottom'); xlim([410 470]); box off;

% Plotting Lever traces for Hits and Misses
figure('Name','Average Lever Traces for Hits & Misses');
subplot(1,2,1);
for i=1:Behaviour.nHit
    plot(Behaviour.hitTrace(i).time-1,Behaviour.hitTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(Behaviour.hitTrace(1).time-parameters.windowBeforePull,mean(horzcat(Behaviour.hitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(10,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in mV)');xlabel('Time (in s)');title('Average Lever Traces for Hits');box off;

subplot(1,2,2);
for i=1:Behaviour.nMiss
    plot(Behaviour.missTrace(i).time-1,Behaviour.missTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(Behaviour.missTrace(1).time-parameters.windowBeforePull,mean(horzcat(Behaviour.missTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(10,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in mV)');xlabel('Time (in s)');title('Average Lever Traces for Misses');box off;

% Plotting Lever traces for Cue Hits and Cue miss 
figure('Name','Average Lever Traces for Cue Hits');
for i=1:Behaviour.nCueHit
    plot(Behaviour.cueHitTrace(i).time-parameters.windowBeforeCue,Behaviour.cueHitTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(Behaviour.cueHitTrace(1).time-parameters.windowBeforeCue,mean(horzcat(Behaviour.cueHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(10,'--.b','Threshold','LabelHorizontalAlignment','left'); 
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(Behaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (in mV)');xlabel('Time (in s)');title('Average Lever Traces for Cue Hits');box off;

%% Reading behaviour data from Intan traces 
IntanBehaviour = readLeverIntan(parameters,LFP.times,Intan.analog_adc_data,Intan.dig_in_data,Behaviour);

% Plotting Lever traces for Hits and Misses
figure('Name','Average Lever Traces for Hits & Misses');
subplot(1,2,1);
for i=1:IntanBehaviour.nHit
    plot(IntanBehaviour.hitTrace(i).time,IntanBehaviour.hitTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.hitTrace(1).time,mean(horzcat(IntanBehaviour.hitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Hits');box off;
subplot(1,2,2);
for i=1:IntanBehaviour.nMiss
    plot(IntanBehaviour.missTrace(i).time,IntanBehaviour.missTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.missTrace(1).time,mean(horzcat(IntanBehaviour.missTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Misses');box off;

% Plotting Lever traces for Cue Hit 
figure('Name','Average Lever Traces for Cue Hits');
for i=1:IntanBehaviour.nCueHit
    plot(IntanBehaviour.cueHitTrace(i).time,IntanBehaviour.cueHitTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.cueHitTrace(1).time,mean(horzcat(IntanBehaviour.cueHitTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
xline(0,'--r','Cue','LabelVerticalAlignment','top');
xline(mean(IntanBehaviour.reactionTime,'all'),'--m','Avg. Reaction Time','LabelVerticalAlignment','top');
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Hits');box off;

% Plotting Lever traces for Cue Miss 
figure('Name','Average Lever Traces for Cue Misses');
for i=1:IntanBehaviour.nCueMiss
    plot(IntanBehaviour.cueMissTrace(i).time,IntanBehaviour.cueMissTrace(i).trace,'Color',[0 0 0 0.2],'LineWidth',1.5);
    hold on;
end
plot(IntanBehaviour.cueMissTrace(1).time,mean(horzcat(IntanBehaviour.cueMissTrace(1:end).trace),2),'Color',[1 0 0 1],'LineWidth',2);
yline(IntanBehaviour.threshold,'--.b','Threshold','LabelHorizontalAlignment','left'); 
xline(0,'--r','Cue','LabelVerticalAlignment','top');
ylabel('Lever deflection (in V)');xlabel('Time (in s)');title('Average Lever Traces for Cue Misses');box off;

%% Generalized Phase 
LFP.xf = bandpass_filter(LFP.LFPdatacube,5,40,4,1000);
[LFP.xgp, LFP.wt] = generalized_phase(LFP.xf,1000,0);
LFP.xfbeta = bandpass_filter(LFP.LFPdatacube,10,30,4,1000);
[LFP.xgpbeta, LFP.wtbeta] = generalized_phase(LFP.xfbeta,1000,0);
LFP.xftheta = bandpass_filter(LFP.LFPdatacube,4,10,4,1000);
[LFP.xgptheta, LFP.wttheta]  = generalized_phase(LFP.xftheta,1000,0);
LFP.xfgamma = bandpass_filter(LFP.LFPdatacube,30,80,4,1000);
[LFP.xgpgamma, LFP.wtgamma]  = generalized_phase(LFP.xfgamma,1000,0);
[parameters.X,parameters.Y] = meshgrid( 1:parameters.cols, 1:parameters.rows );
LFP.xfwide = bandpass_filter(LFP.LFPdatacube,5,90,4,1000);

%% Power Spectrum during task across channels 
[PSDChHit , f] = getAvgPSD(LFP.LFPdatacube,LFP.Fs,IntanBehaviour.cueHitTrace,parameters);
[PSDChMiss , f] = getAvgPSD(LFP.LFPdatacube,LFP.Fs,IntanBehaviour.cueMissTrace,parameters);

avgPSDHit = squeeze(10*log10(mean(PSDChHit,[1 2])));
trialPSDHit = squeeze(10*log10(mean(PSDChHit,2)));
avgPSDMiss = squeeze(10*log10(mean(PSDChMiss,[1 2])));
trialPSDMiss = squeeze(10*log10(mean(PSDChMiss,2)));

figure();
subplot(1,2,1);
plot(f(1:51),trialPSDHit(:,1:51),'Color', [0 0 1 0.1]);
hold on;
plot(f(1:51),avgPSDHit(1:51),'Color', [0 0 1 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Avergage Power Spectral Density for HitTrials');
box off;

subplot(1,2,2);
plot(f(1:51),trialPSDMiss(:,1:51),'Color', [1 0 0 0.1]);
hold on;
plot(f(1:51),avgPSDMiss(1:51),'Color', [1 0 0 1],'LineWidth',1.5);
ylim([0 50]);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Avergage Power Spectral Density for MissTrials');
box off;

trialAvgPSD = 10*log10(squeeze(mean(PSDChHit,1)));
figure();
imagesc(trialAvgPSD(:,1:51)');
colormap("jet");set(gca,'YDir','normal');

%% Wavelet spectrogram
[hitAvgSpectrogram, hitSpectrogramCWT,AvgHitTrace ,fwt] = getAvgSpectogram(LFP.xfwide,LFP.Fs,IntanBehaviour.cueHitTrace,parameters,[5 90]);
[missAvgSpectrogram, missSpectrogramCWT,AvgMissTrace,fwt] = getAvgSpectogram(LFP.xfwide,LFP.Fs,IntanBehaviour.cueMissTrace,parameters,[5 90]);

% Global average spectogram
figure('Name','Trial Averaged Wavelet Spectrogram for Hits & Misses');
subplot(1,2,1);
plotSpectrogram((squeeze(hitAvgSpectrogram)),IntanBehaviour.cueHitTrace(1).time,fwt,'contourf','Wavelet Based Spectrogram for Hits','Time (s)','Frequency (Hz)')
hold on; yyaxis right; box off;
plot(IntanBehaviour.cueHitTrace(1).time,AvgHitTrace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); 
subplot(1,2,2);
plotSpectrogram((squeeze(missAvgSpectrogram)),IntanBehaviour.cueMissTrace(1).time,fwt,'contourf','Wavelet Based Spectrogram for Misses','Time (s)','Frequency (Hz)')
hold on; yyaxis right; box off;
plot(IntanBehaviour.cueMissTrace(1).time,AvgMissTrace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); box off;

% Movement average spectogram 
allAvgSpectogram = mean(cat(1,hitSpectrogramCWT,missSpectrogramCWT));
allAvgLeverTrace = mean(cat(2,horzcat(IntanBehaviour.hitTrace(1:end).trace),horzcat(IntanBehaviour.missTrace(1:end).trace)),2);
figure('Name','Trial Averaged Wavelet Spectrogram for Lever Pull');
plotSpectrogram(squeeze(allAvgSpectogram),IntanBehaviour.hitTrace(1).time-parameters.windowBeforePull,fwt,'Wavelet Based Spectrogram for Lever Pull','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.hitTrace(1).time-parameters.windowBeforePull,allAvgLeverTrace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); box off;

% Plotting for specific trial
trialno = 48;
figure('Name','Spatial Averaged Wavelet Spectrogram for Hits & Misses');
subplot(2,1,1);
plotSpectrogram(squeeze(hitSpectrogramCWT(trialno,:,:)),IntanBehaviour.hitTrace(trialno).time-parameters.windowBeforePull,fwt,'Wavelet based Average Spectogram for Hits','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.hitTrace(trialno).time-parameters.windowBeforePull,IntanBehaviour.hitTrace(trialno).trace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); 
subplot(2,1,2);
trialno = 11;
plotSpectrogram(squeeze(missSpectrogramCWT(trialno,:,:)),IntanBehaviour.missTrace(trialno).time-parameters.windowBeforePull,fwt,'Wavelet based Average Spectogram for Misses','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.missTrace(trialno).time-parameters.windowBeforePull,IntanBehaviour.missTrace(trialno).trace,'-w','LineWidth',2.5);
ylabel('Lever deflection (mV)'); box off;

trialno = 67;
figure('Name','Spatial Averaged Wavelet Spectrogram for Hits & Misses');
subplot(2,1,1);
plotSpectrogram((squeeze(hitSpectrogramCWT(trialno,:,:))),IntanBehaviour.cueHitTrace(trialno).time,fwt,'Wavelet based Average Spectogram for Hits','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.cueHitTrace(trialno).time,squeeze(mean(LFP.xfwide(:,:,IntanBehaviour.cueHitTrace(trialno).LFPIndex),[1 2])),'-w','LineWidth',2);
ylabel('Amplitude (\mu V)'); 
subplot(2,1,2);
trialno = 11;
plotSpectrogram(squeeze(missSpectrogramCWT(trialno,:,:)),IntanBehaviour.missTrace(trialno).time-parameters.windowBeforePull,fwt,'Wavelet based Average Spectogram for Misses','Time (s)','Frequency (Hz)');
hold on; yyaxis right; box off;
plot(IntanBehaviour.missTrace(trialno).time-parameters.windowBeforePull,squeeze(mean(LFP.xfbeta(:,:,IntanBehaviour.missTrace(trialno).LFPIndex),[1 2])),'-w','LineWidth',2);
ylabel('Amplitude (\mu V)');  box off;


%% Initializing plotting options
options.subject = 'W'; % this can be 'W' or 'T' (two marmoset subjects)
options.plot = true; % this option turns plots ON or OFF
options.plot_shuffled_examples = false; % example plots w/channels shuffled in space

%% Wave detection in velocity triggered windows
parameters.spacing = 0.1; % Grid spacing in mm
nShuffle = 1000;
threshold = 99;
trialno = 55;

% Wave detection for wide band
disp('Wave Detection for wide band ...')
rhoThres = getRhoThreshold(LFP.xgp,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
parameters.rhoThres= rhoThres;
Waves.wavesHit = detectWaves(LFP.xf,LFP.xgp,LFP.wt,IntanBehaviour.cueHitTrace,parameters);
if isfield(IntanBehaviour,'cueMissTrace')
    Waves.wavesMiss = detectWaves(LFP.xf,LFP.xgp,LFP.wt,IntanBehaviour.cueMissTrace,parameters);
end


% Wave detection for theta band
disp('Wave Detection for theta band ...')
threshold = 99;
thetarhoThres = getRhoThreshold(LFP.xgptheta,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
parameters.rhoThres = thetarhoThres;
thetaWaves.wavesHit = detectWaves(LFP.xftheta,LFP.xgptheta,LFP.wttheta,IntanBehaviour.cueHitTrace,parameters);
if isfield(IntanBehaviour,'cueMissTrace')
    thetaWaves.wavesMiss = detectWaves(LFP.xftheta,LFP.xgptheta,LFP.wttheta,IntanBehaviour.cueMissTrace,parameters);
end

% Wave detection for beta band
disp('Wave Detection for beta band ...')
threshold = 99;
betarhoThres = getRhoThreshold(LFP.xgpbeta,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
parameters.rhoThres = betarhoThres;
betaWaves.wavesHit= detectWaves(LFP.xfbeta,LFP.xgpbeta,LFP.wtbeta,IntanBehaviour.cueHitTrace,parameters);
if isfield(IntanBehaviour,'cueMissTrace')
    betaWaves.wavesMiss = detectWaves(LFP.xfbeta,LFP.xgpbeta,LFP.wtbeta,IntanBehaviour.cueMissTrace,parameters);
end

% Wave detection for gamma band
disp('Wave Detection for gamma band ...')
threshold = 99;
gammarhoThres = getRhoThreshold(LFP.xgpgamma,IntanBehaviour.cueHitTrace,parameters,nShuffle,trialno,threshold);
parameters.rhoThres = gammarhoThres;    
gammaWaves.wavesHit = detectWaves(LFP.xfgamma,LFP.xgpgamma,LFP.wtgamma,IntanBehaviour.cueHitTrace,parameters);
if isfield(IntanBehaviour,'cueMissTrace')
    gammaWaves.wavesMiss = detectWaves(LFP.xfgamma,LFP.xgpgamma,LFP.wtgamma,IntanBehaviour.cueMissTrace,parameters);
end

%% Wave detecion for entire time 
parameters.rhoThres= rhoThres;
allwaves.LFPIndex = (1:1:size(LFP.LFP,2))';
Wavesall = detectWaves(LFP.xf,LFP.xgp,LFP.wt,allwaves,parameters);
WaveStatsSingle(Wavesall,parameters,0,size(LFP.LFP,2));

%% PLotting to check visually
% trialPlot = 48;
% plot_wave_examples( LFP.xf(:,:,IntanBehaviour.cueHitTrace(trialPlot).LFPIndex(1):IntanBehaviour.cueHitTrace(trialPlot).LFPIndex(end)), ...
%     options, trialPlot, Waves.wavesHit,rhoThres);

% trialPlot = 12;
% plot_wave_examples( LFP.xf(:,:,IntanBehaviour.missTrace(trialPlot).LFPIndex(1):IntanBehaviour.missTrace(trialPlot).LFPIndex(end)), ...
%     options, trialPlot, Waves.wavesMiss,rhoThres);

%% Waves accross trials 
plotOption = 1;
[WaveStats(1)] = getWaveStats(Waves,parameters,plotOption);
[WaveStats(2)] = getWaveStats(thetaWaves,parameters,plotOption);
[WaveStats(3)] = getWaveStats(betaWaves,parameters,plotOption);
[WaveStats(4)] = getWaveStats(gammaWaves,parameters,plotOption);

[WaveStats2(1)] = getInitRewardStats(Waves,parameters,plotOption);
[WaveStats2(2)] = getInitRewardStats(thetaWaves,parameters,plotOption);
[WaveStats2(3)] = getInitRewardStats(betaWaves,parameters,plotOption);
[WaveStats2(4)] = getInitRewardStats(gammaWaves,parameters,plotOption);

%% Percent Phase Locking
[PPLHit] = getPPL(LFP.xgp,IntanBehaviour.cueHitTrace,parameters);
[PPLMiss] = getPPL(LFP.xgp,IntanBehaviour.cueMissTrace,parameters);

figure();
subplot(2,1,1);
imagesc(reshape(PPLHit,[],size(PPLHit,3)));
colormap(hot);
subplot(2,1,2);
imagesc(reshape(PPLMiss,[],size(PPLMiss,3)));
colormap(hot);

figure();
plot(squeeze(nanmean(PPLHit,[1 2])));
hold on;
plot(squeeze(nanmean(PPLMiss,[1 2])));

%% Average PGD accross frequency bands
[PGDfreqHit,PGDfreqMiss] = getPGDFreqBand(LFP.LFPdatacube,IntanBehaviour,1,parameters);

figure();
PGDFreq = [5:5:100];
plot(PGDFreq,PGDfreqHit);
hold on;
plot(PGDFreq,PGDfreqMiss);
xlabel('Frequency');
ylabel('Phase Gradient Directionality');
ylim([0.35 0.7]);
%% Beta event detection 
avgBetaband = mean(LFP.beta_band,1);
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
%% Random code
% %% Plotting video of phase 
% figure(); hold on;
% map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) );axis off;
% x = repmat(1:cols,rows,1); % generate x-coordinates
% y = repmat(1:rows,cols,1)'; % generate y-coordinates
% Ts = 1/1024;
% tvel = ((Ts:Ts:300*Ts))*1000;
% for i=1:numel(goodTrial.xgp)
%     ang = inpaint_nans(angle(goodTrial.xgp(:,:,i)));
% %     ang = inpaint_nans(goodTrial.xf(:,:,i));    
%     pause(0.1);
%     t= num2cell(rad2deg(ang));
%     t = cellfun(@num2str, t, 'UniformOutput',false);
% %     subplot(2,1,1);
%     imagesc(ang); title('Phase map for Electrode Map');
%     text(x(:),y(:),t,'HorizontalAlignment','Center');
% %     ax2 = axes; set( ax2, 'position', [0.02116    0.80976    0.0484    0.1000] ); axis image
% %     [x1,y1] = pol2cart( angle( exp(1i.*linspace(-pi,pi,100)) ), ones( 1, 100 ) );
% %     h3 = cline( x1, y1, linspace(-pi,pi,100) ); axis off; set( h3, 'linewidth', 6 )
% %     t1 = text( 0, 0, 'GP' );
% %     set( t1, 'fontname', 'arial', 'fontsize', 20, 'fontweight', 'bold', 'horizontalalignment', 'center' )
%     
% %     subplot(2,1,2)
% %     plot(Encoder.vel(1+Encoder.velTrig(137)-350:i+Encoder.velTrig(137)-350),'-b','LineWidth',2);xlim([tvel(1) tvel(end)]); ylim([-1 4]);ylabel('Velocity (in cm/s)');xlabel('Time (in ms)');
% end
% 
% %% Plotting video of LFP amplitude
% figure();
% x = repmat(1:cols,rows,1); % generate x-coordinates
% y = repmat(1:rows,cols,1)'; % generate y-coordinates
% tvel = Encoder.timeWindow1;
% for i=1:numel(goodTrial.xgp)
%     amp = inpaint_nans(goodTrial.xf(:,:,i));    
%     pause(0.1);
% %     subplot(2,1,1);
%     imagesc(amp); title('Phase map for Electrode Map');colorbar;
% %     if mod(i,2)==0, continue; end
% %     subplot(2,1,2)
% %     plot(tvel(1:floor(i/2)),Encoder.vel(1,Encoder.trialTime(trialno,1):Encoder.trialTime(trialno,1)+floor(i/2)),'-b','LineWidth',2)
% %     plot(tvel(1:floor(i/2)),'-b','LineWidth',2);
% %     xlim([tvel(1) tvel(end)]); ylim([-1 4]);ylabel('Velocity (in cm/s)');xlabel('Time (in s)');
% end
% %% Finding source points 
% % plot_wave_examples( goodTrial.xf, options, jj, evaluationPoints, source, rho );
% 
% figure();
% for i=1:numel(goodTrial.relTime)
%     quiver(dx(:,:,i),dy(:,:,i));
%     pause(0.1);
%     ylim([0 9]);xlim([0 5]);
% end
% 
% %% 
% 
% dt = 1 / 1024; T = size(LFP.LFP,2) / Fs; time = dt:dt:T;
time = LFP.times(Encoder.trialTime(4,3):Encoder.trialTime(4,4));
xw = squeeze(LFP.xf(4,3,Encoder.trialTime(4,3):Encoder.trialTime(4,4)));
xwRaw = squeeze(LFP.LFPdatacube(4,3,Encoder.trialTime(4,3):Encoder.trialTime(4,4)));
xgp1 = squeeze(LFP.xgp(4,3,Encoder.trialTime(4,3):Encoder.trialTime(4,4)));
% main figure
fg1 = figure; hold on; ax1 = gca; 
plot( time, xwRaw, 'linewidth', 1, 'color', 'k' ); h4 = cline( time, xw, [], angle(xgp1) );
xlim([166 166.5])
set( h4, 'linestyle', '-', 'linewidth', 2  ), axis off
l1 = line( [.1 .2], [-125 -125] ); set( l1, 'linewidth', 4, 'color', 'k' )
l2 = line( [.1 .1], [-125 -75] ); set( l2, 'linewidth', 4, 'color', 'k' )

% inset
map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) )
ax2 = axes; set( ax2, 'position', [0.2116    0.6976    0.0884    0.2000] ); axis image
[x1,y1] = pol2cart( angle( exp(1i.*linspace(-pi,pi,100)) ), ones( 1, 100 ) );
h3 = cline( x1, y1, linspace(-pi,pi,100) ); axis off; set( h3, 'linewidth', 6 )
%%
% text labels
t1 = text( 0, 0, 'GP' );
set( t1, 'fontname', 'arial', 'fontsize', 28, 'fontweight', 'bold', 'horizontalalignment', 'center' )
set( gcf, 'currentaxes', ax1 )
t2 = text( 0.1260, -146.8832, '100 ms' );
set( t2, 'fontname', 'arial', 'fontsize', 24, 'fontweight', 'bold' )
t2 = text( 0.0852, -130.2651, '50 \muV' );
set( t2, 'fontname', 'arial', 'fontsize', 24, 'fontweight', 'bold', 'rotation', 90 )
% 
% %% 
% 
% waveden = zeros(1,size(LFP.LFP,2));
% waveden(evaluationPoints) = 1;
% waveden1 = downsample(waveden,2);
% figure(),plot(waveden1);
% 
% 
% figure();
% plot(Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
% hold on;
% plot(waveden1);
% 
% 
% 
% %meanPhase is the mean phase map over a 5ms window starting at the wave
% %initiaion time
% [rho(jj),~,D,pl] = phase_correlation_distance(meanPhase,source(:,jj), parameters.pixel_spacing);
% %fit the phase (p1) vs distance (D) to a line
% pfit = polyfit(D,pl,1);
% %use phase vs distance to calculate wave speed
% k = abs(pfit(1));%rad/mm, WAVE NUMBER, slope of the lin fit (spatial derivative k=dPhase/dx)
% IF = nanmean(IF,3);%Hz, INSTANTANEOUS FREQUENCY for each electrode in the small time window. f = dPhase/dt
% IF = nanmean(IF(:));%Hz, AVG INSTANTANEOUS FREQUENCY across the grid
% w = IF*2*pi;%(1/s)*(rad) = rad/s, AVG INSTANTANEOUS ANGULAR FREQUENCY across the grid. w = 2pi*f
% speed(jj) = (w/k).*(1/1000); %(rad/s)/(rad/mm) = (mm/s)*(1m/1000mm) = m/s, speed
% 


