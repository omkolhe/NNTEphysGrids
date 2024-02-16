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
%%  PreProcessing
load GridsLowDenNeedle_chanmap.mat;  % load the channel map for the IntanConcatenate function
parameters.rows = 8;  % Number of rows of electrodes on the Grid
parameters.cols = 4;  % Number of colums of electrodes on the Grid
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;
parameters.rewardTrials = 0; % 1 if reward is given , 0 if not
parameters.rewardDistance = 100; % Distance travelled after which reward is given
IntanConcatenate
fpath = Intan.path; % where on disk do you want the analysis? ideally and SSD...

%% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = Ts:Ts:Intan.Tmax;

%% Removing bad channels from impedance values
[Z,Intan.goodChMap,Intan.badChMap] = readImp(electrode_map,5e6);
figure('Name','Impedance Test at 1kHz');boxchart(Z); xlabel('n = ' + string(size(Z,1)));ylabel('Impedance (in \Omega)');set(gca,'xticklabel',{[]})
%Intan = removeBadCh(Intan,Intan.badCh);

%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(Intan.allIntan,10000);
% LFP = bestLFP(LFP);
% LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
% LFPplot(LFP);
LFP = createDataCube(LFP,parameters.rows,parameters.cols,Intan.badChMap); % Creating datacube

%% Spikes 
Spikes.hpSignal = bandpass_filter_matrix(double(Intan.allIntan),300,3000,10,10000);
Spikes.hpSignal(Intan.badChMap,:) = NaN;
Spikes = CAR(Spikes);

figure,stack_plot(Spikes.hpSpikes(:,20000:30000),1,4,10000);
figure,stack_plot(Spikes.hpSignal(:,1:10000),1,4,10000);


%% Loading Encoder Data
[Encoder] = readPos(parameters);
figure('Name','Velocity');plot(Encoder.time,Encoder.vel,'LineWidth',1.5);ylim([-20 20]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
Encoder.vel = abs(Encoder.vel);
%Encoder.vel(Encoder.vel<0) = 0;

%% Find velocity triggers before velocity reaches 2cm/s
Encoder = detectVelTrigStart(Encoder,2,0.05,50,100,LFP.times); % For calculating initiation
Encoder = detectVelTrigStop(Encoder,2,0.05,50,100,LFP.times); % For calculating stopping 
figure('Name','Velocity');plot(Encoder.time,Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
hold on;xline(Encoder.velTrig(1,:)/Encoder.fs,'b');xline(Encoder.velTrigStop(1,:)/Encoder.fs,'r');


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

%% Segementing trial windows
windowBeforeTrig = 2; % in seconds
windowAfterTrig = 2; % in seconds
Encoder = getVelTrigTrials(Encoder,windowBeforeTrig,windowAfterTrig,LFP.Fs); % Selecting good trials for initiation
Encoder = getVelTrigTrialsStop(Encoder,windowBeforeTrig,windowAfterTrig,LFP.Fs); % Selecting good trials for stopping
Encoder.timeWindow1 = -1*windowBeforeTrig:1/Encoder.fs:windowAfterTrig-1/Encoder.fs; % Time series for plotting
Encoder.timeWindow2 = -1*windowBeforeTrig:1/LFP.Fs:windowAfterTrig-1/LFP.Fs; % Time series for plotting

figure('Name','Velocity Triggers trials');
subplot(1,2,1);
for i=1:Encoder.nTrig
    a = plot(Encoder.timeWindow1,Encoder.vel(Encoder.trialTime(i,1):Encoder.trialTime(i,2)),'LineWidth',1.5);
    Encoder.velTrial(i,:) = Encoder.vel(Encoder.trialTime(i,1):Encoder.trialTime(i,2));
    ylim([0 20]);xlabel('Time (in s)');ylabel('Velocity in cm/s');%yline([2 -2]);
    hold on;xline(0,'--');title("Velocity for Initiation Trials");box off;
    text(0.75, 19, 'n = ' + string(Encoder.nTrig),'HorizontalAlignment', 'right','VerticalAlignment', 'middle','Color', 'black');
end
subplot(1,2,2);
for i=1:Encoder.nTrigStop
    plot(Encoder.timeWindow1,Encoder.vel(Encoder.trialTimeStop(i,1):Encoder.trialTimeStop(i,2)),'LineWidth',1.5);
    Encoder.velTrialStop(i,:) = Encoder.vel(Encoder.trialTimeStop(i,1):Encoder.trialTimeStop(i,2));
    ylim([0 20]);xlabel('Time (in s)');ylabel('Velocity in cm/s');%yline([2 -2]);
    hold on;xline(0,'--');title("Velocity for Termination Trials"); box off;
    text(0.75, 19, 'n = ' + string(Encoder.nTrigStop),'HorizontalAlignment', 'right','VerticalAlignment', 'middle','Color', 'black');
end

%% Segmenting Encoder data into Motion States
Encoder = getMotionStates(Encoder,LFP.times);
figure('Name','Motion states');
yyaxis left;plot(Encoder.time,Encoder.discVel,'LineWidth',1.5);ylabel('Velocity (in cm/s)');xlabel('Time (in s)'); hold on
yyaxis right; plot(Encoder.time,Encoder.state,'-r','LineWidth',1.5);ylabel('State'); 
dim = [0.75 0.5 0.25 0.3];
str = {'0 - Rest','1 - Initiation','2 - Run','3 - Termination'};
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%% Wavelet spectrogram
[globalAvgSpectrogram, avgSpectrogramCWT,avgVel,fwt] = getAvgSpectogram(LFP.LFPdatacube,LFP,Encoder,Encoder.trialTime(:,:),Encoder.velTrial,parameters,[5 90]);
[globalAvgSpectrogramStop, avgSpectrogramCWTStop, avgVelStop,fwtStop] = getAvgSpectogram(LFP.LFPdatacube,LFP,Encoder,Encoder.trialTimeStop(:,:),Encoder.velTrialStop,parameters,[5 45]);

trialno = 2;
figure('Name','Trial Averaged Wavelet Spectrogram for Motion Inititation');
globalAvgVel = interp(mean(Encoder.velTrial(trialno,:),1),LFP.Fs/Encoder.fs);
globalAvgLFP = squeeze(mean(LFP.xfbeta(:,:,Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4)),[1 2]));
imagesc(Encoder.timeWindow2,fwt,squeeze(globalAvgSpectrogram));colormap('jet');set(gca,'YDir','normal');title('Wavelet based Average Spectogram');ylabel('Frequency (Hz)');xlabel('Time (s)');
c=colorbar;ylabel(c, 'Relative Power to white noise','FontSize',10);hold on; yyaxis right; box off;
plot(Encoder.timeWindow2,globalAvgLFP,'-r','LineWidth',2.5);
% plot(Encoder.timeWindow2,globalAvgVel,'-w','LineWidth',2.5);ylabel('Velocity (cm/s)');xlabel('Time (ms)');

figure('Name','Trial Averaged Wavelet Spectrogram for Motion Termination');
globalAvgVelStop = interp(mean(Encoder.velTrialStop(1,:),1),LFP.Fs/Encoder.fs);
imagesc(Encoder.timeWindow2,fwtStop,squeeze(globalAvgSpectrogramStop));colormap('jet');set(gca,'YDir','normal');title('Wavelet based Average Spectogram');ylabel('Frequency (Hz)');xlabel('Time (s)');
c=colorbar;ylabel(c, 'Relative Power to white noise','FontSize',10);hold on;yyaxis right; box off;
plot(Encoder.timeWindow2,globalAvgVelStop,'-w','LineWidth',2.5);ylabel('Velocity (cm/s)');xlabel('Time (ms)');

%% Initializing plotting options
options.subject = 'W'; % this can be 'W' or 'T' (two marmoset subjects)
options.plot = true; % this option turns plots ON or OFF
options.plot_shuffled_examples = false; % example plots w/channels shuffled in space

%% Wave detection in velocity triggered windows
parameters.spacing = 0.1; % Grid spacing in mm
nShuffle = 10000;
threshold = 99;
trialno = 1;

% Wave detection for wide band
disp('Wave Detection for wide band ...')
% rhoThres = getRhoThreshold(LFP.xgp,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres= rhoThres;
Waves.wavesStart = detectWaves(LFP.xf,LFP.xgp,LFP.wt,Encoder.trialTime,parameters);
Waves.wavesStop = detectWaves(LFP.xf,LFP.xgp,LFP.wt,Encoder.trialTimeStop,parameters);

% Wave detection for theta band
disp('Wave Detection for theta band ...')
threshold = 99;
% thetarhoThres = getRhoThreshold(LFP.xgptheta,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres = thetarhoThres;
thetaWaves.wavesStart = detectWaves(LFP.xftheta,LFP.xgptheta,LFP.wttheta,Encoder.trialTime,parameters);
thetaWaves.wavesStop = detectWaves(LFP.xftheta,LFP.xgptheta,LFP.wttheta,Encoder.trialTimeStop,parameters);

% Wave detection for beta band
disp('Wave Detection for beta band ...')
threshold = 99.9;
% betarhoThres = getRhoThreshold(LFP.xgpbeta,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres = betarhoThres;
betaWaves.wavesStart = detectWaves(LFP.xfbeta,LFP.xgpbeta,LFP.wtbeta,Encoder.trialTime,parameters);
betaWaves.wavesStop = detectWaves(LFP.xfbeta,LFP.xgpbeta,LFP.wtbeta,Encoder.trialTimeStop,parameters);

% Wave detection for gamma band
disp('Wave Detection for gamma band ...')
threshold = 99.9;
% gammarhoThres = getRhoThreshold(LFP.xgpgamma,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres = gammarhoThres;    
gammaWaves.wavesStart = detectWaves(LFP.xfgamma,LFP.xgpgamma,LFP.wtgamma,Encoder.trialTime,parameters);
gammaWaves.wavesStop = detectWaves(LFP.xfgamma,LFP.xgpgamma,LFP.wtgamma,Encoder.trialTimeStop,parameters);

%% PLotting to check visually
trialPlot = 1;
plot_wave_examples( LFP.xf(:,:,Encoder.trialTime(trialPlot,3):Encoder.trialTime(trialPlot,4)), ...
    options, trialPlot, Waves.wavesStart,rhoThres);

plot_wave_examples( LFP.xf(:,:,Encoder.trialTimeStop(trialPlot,3):Encoder.trialTimeStop(trialPlot,4)), ...
    options, trialPlot, Waves.wavesStop,rhoThres);

plot_wave_examples( LFP.xfbeta(:,:,Encoder.trialTime(trialPlot,3):Encoder.trialTime(trialPlot,4)), ...
    options, trialPlot, betaWaves,betarhoThres);

%% Waves accross trials 
plotOption = 1;
[WaveStats(1)] = getWaveStats(Waves,parameters,plotOption);
[WaveStats(2)] = getWaveStats(thetaWaves,parameters,plotOption);
[WaveStats(3)] = getWaveStats(betaWaves,parameters,plotOption);
[WaveStats(4)] = getWaveStats(gammaWaves,parameters,plotOption);

%% Wave detecion for entire time 
parameters.rhoThres= rhoThres;
Wavesall = detectWaves(LFP.xgp,LFP.wt,[0,0,1,size(LFP.xgp,3)],parameters);
WaveStatsSingle(Wavesall,parameters,1);

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

%% Wave detection for Rest, Run, Init and Term

% Wave detection for wide band
disp('Wave Detection for wide band ...')
% rhoThres = getRhoThreshold(LFP.xgp,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres= rhoThres;
Waves.wavesRest = detectWaves(LFP.xf,LFP.xgp,LFP.wt,Encoder.restTrialTime,parameters);
Waves.wavesRun = detectWaves(LFP.xf,LFP.xgp,LFP.wt,Encoder.runTrialTime,parameters);
Waves.wavesInit = detectWaves(LFP.xf,LFP.xgp,LFP.wt,Encoder.initTrialTime,parameters);
Waves.wavesTerm = detectWaves(LFP.xf,LFP.xgp,LFP.wt,Encoder.termTrialTime,parameters);

% Wave detection for wide band
disp('Wave Detection for theta band ...')
parameters.rhoThres= thetarhoThres;
thetaWaves.wavesRest = detectWaves(LFP.xftheta,LFP.xgptheta,LFP.wttheta,Encoder.restTrialTime,parameters);
thetaWaves.wavesRun = detectWaves(LFP.xftheta,LFP.xgptheta,LFP.wttheta,Encoder.runTrialTime,parameters);
thetaWaves.wavesInit = detectWaves(LFP.xftheta,LFP.xgptheta,LFP.wttheta,Encoder.initTrialTime,parameters);
thetaWaves.wavesTerm = detectWaves(LFP.xftheta,LFP.xgptheta,LFP.wttheta,Encoder.termTrialTime,parameters);

% Wave detection for beta band
disp('Wave Detection for beta band ...')
parameters.rhoThres= betarhoThres;
betaWaves.wavesRest = detectWaves(LFP.xfbeta,LFP.xgpbeta,LFP.wtbeta,Encoder.restTrialTime,parameters);
betaWaves.wavesRun = detectWaves(LFP.xfbeta,LFP.xgpbeta,LFP.wtbeta,Encoder.runTrialTime,parameters);
betaWaves.wavesInit = detectWaves(LFP.xfbeta,LFP.xgpbeta,LFP.wtbeta,Encoder.initTrialTime,parameters);
betaWaves.wavesTerm = detectWaves(LFP.xfbeta,LFP.xgpbeta,LFP.wtbeta,Encoder.termTrialTime,parameters);

% Wave detection for gamma band
disp('Wave Detection for gamma band ...')
% rhoThres = getRhoThreshold(LFP.xgp,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres= gammarhoThres;
gammaWaves.wavesRest = detectWaves(LFP.xfgamma,LFP.xgpgamma,LFP.wtgamma,Encoder.restTrialTime,parameters);
gammaWaves.wavesRun = detectWaves(LFP.xfgamma,LFP.xgpgamma,LFP.wtgamma,Encoder.runTrialTime,parameters);
gammaWaves.wavesInit = detectWaves(LFP.xfgamma,LFP.xgpgamma,LFP.wtgamma,Encoder.initTrialTime,parameters);
gammaWaves.wavesTerm = detectWaves(LFP.xfgamma,LFP.xgpgamma,LFP.wtgamma,Encoder.termTrialTime,parameters);


% Getting stats 
plotOption = 0;
[WaveStatsStates(1)] = getWaveStatsStates(Waves,parameters,plotOption);
[WaveStatsStates(2)] = getWaveStatsStates(thetaWaves,parameters,plotOption);
[WaveStatsStates(3)] = getWaveStatsStates(betaWaves,parameters,plotOption);
[WaveStatsStates(4)] = getWaveStatsStates(gammaWaves,parameters,plotOption);
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
%%

% speedComb = [WaveStatsT4(1).speed WaveStatsT3(1).speed WaveStatsT2(1).speed];
% avgSpeed = mean(speedComb);
% 
% speedCombStop = [WaveStatsT4(1).speedStop WaveStatsT3(1).speedStop WaveStatsT2(1).speedStop];
% avgSpeedStop = mean(speedCombStop);

speedComb = [WaveStatsT4(3).speed WaveStatsT3(3).speed WaveStatsT2(3).speed];
avgSpeed = mean(speedComb);

speedCombStop = [WaveStatsT4(4).speed WaveStatsT3(4).speed WaveStatsT2(4).speed];
avgSpeedStop = mean(speedCombStop);

% Perform the t-test.
[p, t] = ranksum(speedComb, speedCombStop);
% Print the results.
disp('Wave Speed')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);

figure('Name','Wave speeds in motion Initiation and termination');
group = [ones(size(speedComb')); 2.*ones(size(speedCombStop'))]; 
boxplot([speedComb';speedCombStop'],group,'BoxStyle','filled','PlotStyle','compact');
set(gca,'xtick',1:2,'XTickLabel',{'Initiation','Termination'});box off;
ylabel('Wave speed in cm/s');


% dirComb = [WaveStatsT2(1).velDir WaveStatsT3(1).velDir]; %WaveStatsT3(1).speed WaveStatsT2(1).speed];
% avgDir = mean(dirComb);
% 
% dirCombStop = [WaveStatsT2(1).velDirStop WaveStatsT3(1).velDirStop]; %WaveStatsT3(1).speedStop WaveStatsT2(1).speedStop];
% avgvelDirStop = mean(dirCombStop);

dirComb = [ WaveStatsT4(3).velDir]; %WaveStatsT3(1).speed WaveStatsT2(1).speed];
avgDir = mean(dirComb);

dirCombStop = [ WaveStatsT4(4).velDir]; %WaveStatsT3(1).speedStop WaveStatsT2(1).speedStop];
avgvelDirStop = mean(dirCombStop);


[p, t] = ranksum(dirComb, dirCombStop);
% Print the results.
disp('Wave Direction')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


figure('Name','Polar Histogram for wave direction in Motion Initiation and Termination');
subplot(2,1,1);
polarhistogram(dirComb,30); box off;
title('Wave Direction : Motion Initiation');
subplot(2,1,2);
polarhistogram(dirCombStop,30); box off;
title('Wave Direction : Motion Termination');

figure('Name','Wave Direction in motion Initiation and termination');
group = [ones(size(dirComb'));2.*ones(size(dirCombStop'))];
boxplot([rad2deg(dirComb)';rad2deg(dirCombStop)'],group,'BoxStyle','filled','PlotStyle','compact');
set(gca,'XTickLabel',{'Initiation','Termination'});


speedCombRest = [WaveStatsStatesT4(3).speedRest WaveStatsStatesT3(3).speedRest WaveStatsStatesT2(3).speedRest];
avgSpeedRest = mean(speedCombRest);
speedCombRun = [WaveStatsStatesT4(3).speedRun WaveStatsStatesT3(3).speedRun WaveStatsStatesT2(3).speedRun];
avgSpeedRun = mean(speedCombRun);
speedCombInit = [WaveStatsStatesT4(3).speedInit WaveStatsStatesT3(3).speedInit WaveStatsStatesT2(3).speedInit];
avgSpeedInit = mean(speedCombInit);
speedCombTerm = [WaveStatsStatesT4(3).speedTerm WaveStatsStatesT3(3).speedTerm WaveStatsStatesT2(3).speedTerm];
avgSpeedTerm = mean(speedCombTerm);

% Perform the t-test.
[p, t] = ranksum(speedCombRun , speedCombTerm);
% Print the results.
disp('Wave Speed')
disp('h-statistic:');
disp(t);
disp('p-value:');
disp(p);


figure('Name','Histogram of wave speeds');
subplot(4,1,1);
histfit(speedCombRest,100,'kernel');
xline(avgSpeedRest,'-r',{'Mean speed = ' num2str(avgSpeedRest) ' cm/s'});box off;
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Rest');xlim([0 inf]);
subplot(4,1,2);
histfit(speedCombInit,100,'kernel');
xline(avgSpeedInit,'-r',{'Mean speed = ' num2str(avgSpeedInit) ' cm/s'});box off;
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Initiation');xlim([0 inf]);
subplot(4,1,3);
histfit(speedCombRun,100,'kernel');
xline(avgSpeedRun,'-r',{'Mean speed = ' num2str(avgSpeedRun) ' cm/s'});box off;
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Run');xlim([0 inf]);
subplot(4,1,4);
histfit(speedCombTerm,100,'kernel');
xline(avgSpeedTerm,'-r',{'Mean speed = ' num2str(avgSpeedTerm) ' cm/s'});box off;
xlabel('Wave speed in cm/s');ylabel('Frequency');title('Wave Speed : Motion Term');xlim([0 inf]);

figure('Name','Wave speeds in Rest, Initiation, Running and Termination');
group = [ones(size(speedCombRest')); 2.*ones(size(speedCombInit')); 3.*ones(size(speedCombRun')); 4.*ones(size(speedCombTerm'))];
boxplot([speedCombRest';speedCombInit';speedCombRun';speedCombTerm'],group,'BoxStyle','filled','PlotStyle','compact');
set(gca,'xtick',1:4,'XTickLabel',{'Rest','Initiation','Run','Termination'});box off;
ylabel('Wave speed in cm/s');



