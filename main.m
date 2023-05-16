% clear; clc; 
% close all;
format compact;
set(0,'DefaultFigureWindowStyle','normal')
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
parameters.cols = 4;  % Numever of colums of electrodes on the Grid
IntanConcatenate
fpath = Intan.path; % where on disk do you want the analysis? ideally and SSD...

%% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = Ts:Ts:Intan.Tmax;

%% Removing bad channels from impedance values
[Z,Intan.goodChMap,Intan.badChMap] = readImp(electrode_map,100e6);
figure('Name','Impedance Test at 1kHz');boxchart(Z); xlabel('n = ' + string(size(Z,1)));ylabel('Impedance (in \Omega)');set(gca,'xticklabel',{[]})
%Intan = removeBadCh(Intan,Intan.badCh);

%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(Intan.allIntan,8000);
LFP = bestLFP(LFP);
%LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
LFPplot(LFP);
LFP = createDataCube(LFP,parameters.rows,parameters.cols,Intan.goodChMap); % Creating datacube

%% Loading Encoder Data
[Encoder.pos, Encoder.vel, Encoder.time, Encoder.fs] = readPos();
figure('Name','Velocity');plot(Encoder.time,Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
Encoder.vel = abs(Encoder.vel);
%Encoder.vel(Encoder.vel<0) = 0;

%% Find velocity triggers before velocity reaches 2cm/s
Encoder = detectVelTrig(Encoder,2,0.05,50,100,LFP.times);
figure('Name','Velocity');plot(Encoder.time,Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
hold on;xline(Encoder.velTrig(1,:)/Encoder.fs);

%% Generalized Phase 
LFP.xf = bandpass_filter(LFP.LFPdatacube,5,40,4,1000);
[LFP.xgp, LFP.wt] = generalized_phase(LFP.xf,1000,0);
LFP.xfbeta = bandpass_filter(LFP.LFPdatacube,10,30,4,1000);
[LFP.xgpbeta, LFP.wtbeta] = generalized_phase(LFP.xfbeta,1000,0);
LFP.xftheta = bandpass_filter(LFP.LFPdatacube,4,10,4,1000);
[LFP.xgptheta, LFP.wttheta]  = generalized_phase(LFP.xftheta,1000,0);
LFP.xfgamma = bandpass_filter(LFP.LFPdatacube,30,80,4,1000);
[LFP.xgpgamma, LFP.wtgamma]  = generalized_phase(LFP.xfgamma,1000,0);
[X,Y] = meshgrid( 1:parameters.cols, 1:parameters.rows );

%% Segementing trial windows
windowBeforeTrig = 0.6; % in seconds
windowAfterTrig = 0.4; % in seconds
Encoder = getVelTrigTrials(Encoder,windowBeforeTrig,windowAfterTrig,LFP.Fs);

Encoder.timeWindow1 = -1*windowBeforeTrig:1/Encoder.fs:windowAfterTrig-1/Encoder.fs; % Time series for plotting
Encoder.timeWindow2 = -1*windowBeforeTrig:1/LFP.Fs:windowAfterTrig-1/LFP.Fs; % Time series for plotting
figure('Name','Velocity Triggers trials');
for i=1:Encoder.nTrig
    plot(Encoder.timeWindow1,Encoder.vel(Encoder.trialTime(i,1):Encoder.trialTime(i,2)),'LineWidth',1.5)
    Encoder.velTrial(i,:) = Encoder.vel(Encoder.trialTime(i,1):Encoder.trialTime(i,2));
    ylim([0 20]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
    hold on;xline(0);
end

%% Initializing plotting options
options.subject = 'W'; % this can be 'W' or 'T' (two marmoset subjects)
options.plot = true; % this option turns plots ON or OFF
options.plot_shuffled_examples = false; % example plots w/channels shuffled in space

%% Wave detection in velocity triggered windows
parameters.spacing = 0.1; % Grid spacing in mm
parameters.X = X;
parameters.Y = Y;
nShuffle = 10000;
threshold = 99;
trialno = 1;

% Wave detection for wide band
disp('Wave Detection for wide band ...')
rhoThres = getRhoThreshold(LFP.xgp,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres= rhoThres;
Waves = detectWaves(LFP.xgp,LFP.wt,Encoder,parameters);

% Wave detection for theta band
disp('Wave Detection for theta band ...')
threshold = 99;
thetarhoThres = getRhoThreshold(LFP.xgptheta,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres = thetarhoThres;
thetaWaves = detectWaves(LFP.xgptheta,LFP.wttheta,Encoder,parameters);

% Wave detection for beta band
disp('Wave Detection for beta band ...')
threshold = 99.9;
betarhoThres = getRhoThreshold(LFP.xgpbeta,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres = betarhoThres;
betaWaves = detectWaves(LFP.xgpbeta,LFP.wtbeta,Encoder,parameters);

% Wave detection for gamma band
disp('Wave Detection for gamma band ...')
threshold = 99.9;
gammarhoThres = getRhoThreshold(LFP.xgpgamma,Encoder,parameters,nShuffle,trialno,threshold);
parameters.rhoThres = gammarhoThres;    
gammaWaves = detectWaves(LFP.xgpgamma,LFP.wtgamma,Encoder,parameters);
%% PLotting to check visually
trialPlot = 5;
plot_wave_examples( LFP.xf(:,:,Encoder.trialTime(trialPlot,3):Encoder.trialTime(trialPlot,4)), ...
    options, trialPlot, Waves,rhoThres);

plot_wave_examples( LFP.xfbeta(:,:,Encoder.trialTime(trialPlot,3):Encoder.trialTime(trialPlot,4)), ...
    options, trialPlot, betaWaves,betarhoThres);

%% Waves accross trials 

[WaveStats(1)] = getWaveStats(Waves,parameters,1);
[WaveStats(2)] = getWaveStats(thetaWaves,parameters,0);
[WaveStats(3)] = getWaveStats(betaWaves,parameters,1);
[WaveStats(4)] = getWaveStats(gammaWaves,parameters,0);

%% Wavelet spectrogram
for trialno = 1:size(Encoder.velTrig,2)
    goodTrial.xf = LFP.xf(:,:,Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4));
    goodTrial.xgp = LFP.xgp(:,:,Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4));
    goodTrial.time = LFP.times(Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4));
    goodTrial.relTime = Encoder.timeWindow2;
    figure();
    for i=1:32
        subplot(parameters.rows,parameters.cols,i);
        if ismember(i,Intan.badChMap), continue; end
        calCWTSpectogram(squeeze(goodTrial.xf(floor((i-1)/parameters.cols)+1,mod(i-1,parameters.cols)+1,:)),goodTrial.relTime,1024,20,[1 45],0);
    %     [s,f,t,ps,fc,tc] = spectrogram(squeeze(goodTrial.xf(floor((i-1)/cols)+1,mod(i-1,cols)+1,:)),25,16,1024,1024,'yaxis','onesided');
    %     imagesc(t,f(10:45),abs(ps(10:45,:)));hold on; xline(301/1024,'-r', 'LineWidth',2);set(gca,'YDir','normal') 
    end
    
    % Average spectrogram across all channels
    for i=1:parameters.rows
        for j=1:parameters.cols
            [spectrogramCh((i-1)*parameters.cols + j,:,:) ,fwt] = calCWTSpectogram(squeeze(goodTrial.xf(parameters.rows,parameters.cols,:)),goodTrial.relTime,1024,20,[1 45],0);
        end
    end
    avgSpectrogramCWT(trialno,:,:) = mean(spectrogramCh,1);
end
figure();
globalAvgSpeectrogram = mean(avgSpectrogramCWT,1);
globalAvgVel = mean(Encoder.velTrial,1);
subplot(2,1,1);
imagesc(goodTrial.relTime,fwt,squeeze(globalAvgSpeectrogram));colormap('jet');set(gca,'YDir','normal');title('Wavelet based Average Spectogram');ylabel('Frequency (Hz)');xlabel('Time (s)');
c=colorbar;ylabel(c, 'Relative Power to white noise','FontSize',10);
subplot(2,1,2);
plot(Encoder.timeWindow1,globalAvgVel,'-r','LineWidth',2);
ylabel('Velocity (cm/s)');xlabel('Time (ms)');

%     subplot(2,1,1);
%     imagesc(goodTrial.relTime,fwt,squeeze(avgSpectrogramCWT(1,:,:)));colormap('jet');set(gca,'YDir','normal');title('Wavelet based Average Spectogram');ylabel('Frequency (Hz)');xlabel('Time (s)');
%     c=colorbar;ylabel(c, 'Relative Power to white noise','FontSize',10);
%     subplot(2,1,2);
%     plot(Encoder.timeWindow1,Encoder.vel(Encoder.trialTime(trialno,1):Encoder.trialTime(trialno,2)),'-r','LineWidth',2);
%     ylabel('Velocity (cm/s)');xlabel('Time (ms)');

%% PLotting LFP 
figure();
for i=1:32
    subplot(parameters.rows,parameters.cols,i);
    if ismember(i,Intan.badChMap), continue; end
    plot(goodTrial.relTime,squeeze(goodTrial.xf(floor((i-1)/parameters.cols)+1,mod(i-1,parameters.cols)+1,:))');xline(0,'-r');
end
%% Calculating waves without windowing
p = LFP.xgp(:,:,:);
wt = LFP.wt(:,:,:);
%p = arrayfun(@(jj) inpaint_nans(p(:,:,jj)),1:size(p,3));
evaluationPoints = find_evaluation_points(p,pi,0.2);
%plot_evaluation_points( p, evaluationPoints );
[pm,pd,dx,dy] = phase_gradient_complex_multiplication( p, spacing);
% Phase gradient directionality 
PGD = phase_gradient_directionality(pm,dx,dy);
% Wavelength
wl = 1./abs(pm);
% Instantaneous speed   
s = instantaneous_speed(wt,pm);
% divergence calculation
source = find_source_points(evaluationPoints, X, Y, dx, dy );
% phase correlation with distance (\rho_{\phi,d} measure)
rho = zeros( 1, length(evaluationPoints) );
for jj = 1:length(evaluationPoints)
    ph = angle( p(:,:,evaluationPoints(jj)) );
    rho(jj) = phase_correlation_distance( ph,source(:,jj), spacing );
end
indxBadWave = find(rho<rhoThres);
evaluationPoints(indxBadWave) = [];
source(:,indxBadWave) = [];
rho(indxBadWave) = [];
nWaves = size(evaluationPoints,2);

figure('Name','Velocity');plot(Encoder.time,Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
hold on;xline(evaluationPoints/LFP.Fs);

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
% % dt = 1 / 1024; T = size(LFP.LFP,2) / Fs; time = dt:dt:T;
% time = LFP.times(500:6000);
% xw = squeeze(LFP.xf(4,3,500:6000));
% xwRaw = squeeze(LFP.LFPdatacube(4,3,500:6000));
% xgp1 = squeeze(LFP.xgp(4,3,500:6000));
% % main figure
% fg1 = figure; hold on; ax1 = gca; 
% plot( time, xw, 'linewidth', 2, 'color', 'k' ); h4 = cline( time, xw, [], angle(xgp1) );
% xlim([1 2])
% set( h4, 'linestyle', '-', 'linewidth', 2  ), axis off
% l1 = line( [.1 .2], [-125 -125] ); set( l1, 'linewidth', 4, 'color', 'k' )
% l2 = line( [.1 .1], [-125 -75] ); set( l2, 'linewidth', 4, 'color', 'k' )
% 
% % inset
% map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) )
% ax2 = axes; set( ax2, 'position', [0.2116    0.6976    0.0884    0.2000] ); axis image
% [x1,y1] = pol2cart( angle( exp(1i.*linspace(-pi,pi,100)) ), ones( 1, 100 ) );
% h3 = cline( x1, y1, linspace(-pi,pi,100) ); axis off; set( h3, 'linewidth', 6 )
% %%
% % text labels
% t1 = text( 0, 0, 'GP' );
% set( t1, 'fontname', 'arial', 'fontsize', 28, 'fontweight', 'bold', 'horizontalalignment', 'center' )
% set( gcf, 'currentaxes', ax1 )
% t2 = text( 0.1260, -146.8832, '100 ms' );
% set( t2, 'fontname', 'arial', 'fontsize', 24, 'fontweight', 'bold' )
% t2 = text( 0.0852, -130.2651, '50 \muV' );
% set( t2, 'fontname', 'arial', 'fontsize', 24, 'fontweight', 'bold', 'rotation', 90 )
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
