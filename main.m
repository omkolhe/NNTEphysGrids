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
spacing = 0.1; % Grid spacing in mm
rows = 8;  % Number of rows of electrodes on the Grid
cols = 4;  % Numevbr of colums of electrodes on the Grid
IntanConcatenate
fpath    = Intan.path; % where on disk do you want the analysis? ideally and SSD...

%% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = Ts:Ts:Intan.Tmax;

%% Removing bad channels from impedance values
[Z,Intan.goodChMap,Intan.badChMap] = readImp(electrode_map,100e6);
figure('Name','Impedance Test at 1kHz');boxchart(Z); xlabel('n = ' + string(size(Z,1)));ylabel('Impedance (in \Omega)');set(gca,'xticklabel',{[]})
% Intan = removeBadCh(Intan,Intan.badCh);

%% LFP
set(0,'DefaultFigureWindowStyle','normal')
% LFP = fastpreprocess_filtering(flip(Intan.allIntan,1),8192) ; %Only run for PFF data
LFP = fastpreprocess_filtering(Intan.allIntan,8000);
LFP = bestLFP(LFP);
LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
LFPplot(LFP);
LFP = createDataCube(LFP,rows,cols,Intan.goodChMap); % Creating datacube

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
[X,Y] = meshgrid( 1:cols, 1:rows );

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

%% Wavelet spectrogram
for trialno = 1:size(Encoder.velTrig,2)
    goodTrial.xf = LFP.xf(:,:,Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4));
    goodTrial.xgp = LFP.xgp(:,:,Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4));
    goodTrial.time = LFP.times(Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4));
    goodTrial.relTime = Encoder.timeWindow2;
    figure();
    for i=1:32
        subplot(rows,cols,i);
        if ismember(i,Intan.badChMap), continue; end
        calCWTSpectogram(squeeze(goodTrial.xf(floor((i-1)/cols)+1,mod(i-1,cols)+1,:)),goodTrial.relTime,1024,20,[1 45],0);
    %     [s,f,t,ps,fc,tc] = spectrogram(squeeze(goodTrial.xf(floor((i-1)/cols)+1,mod(i-1,cols)+1,:)),25,16,1024,1024,'yaxis','onesided');
    %     imagesc(t,f(10:45),abs(ps(10:45,:)));hold on; xline(301/1024,'-r', 'LineWidth',2);set(gca,'YDir','normal') 
    end
    
    % Average spectrogram across all channels
    for i=1:rows
        for j=1:cols
            [spectrogramCh((i-1)*cols + j,:,:) ,fwt] = calCWTSpectogram(squeeze(goodTrial.xf(rows,cols,:)),goodTrial.relTime,1024,20,[1 45],0);
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
    subplot(rows,cols,i);
    if ismember(i,Intan.badChMap), continue; end
    plot(goodTrial.relTime,squeeze(goodTrial.xf(floor((i-1)/cols)+1,mod(i-1,cols)+1,:))');xline(0,'-r');
end

%% rho for shuffled data
nShuffle = 10000;
for ii=1:1%size(Encoder.velTrig,2)
    p = LFP.xgp(:,:,Encoder.trialTime(ii,3):Encoder.trialTime(ii,4));
    evaluationPoints = find_evaluation_points(p,pi,spacing);
    rho = zeros( nShuffle, length(evaluationPoints) );
    for kk=1:nShuffle
        if kk==1
            pShuffle = p;
        else
            pShuffle = shuffle_channels(p);
        end
        [pm,pd,dx,dy] = phase_gradient_complex_multiplication( pShuffle, spacing );
        source = find_source_points( evaluationPoints, X, Y, dx, dy );
        for jj = 1:length(evaluationPoints)
            ph = angle( pShuffle(:,:,evaluationPoints(jj)) );
            rho(kk,jj) = phase_correlation_distance( ph, source(:,jj), spacing );
        end
    end
end

% Plotting rho distribution
nEvalPoints = size(rho,2);
rho1 = reshape(rho',[1,nShuffle*nEvalPoints]);
rhoThres = prctile(rho1,99.9);
figure();histogram(rho1(1:nEvalPoints),'FaceColor','r'); hold on;
histogram(rho1(nEvalPoints+1:end),'FaceColor','b'); 
xline(rhoThres,'-r',{'99.9 Percentile'});
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
%% Pre allocating and initializing plotting options
options.subject = 'W'; % this can be 'W' or 'T' (two marmoset subjects)
options.plot = true; % this option turns plots ON or OFF
options.plot_shuffled_examples = false; % example plots w/channels shuffled in space
% Waves = struct(size(Encoder.velTrig,2),1);

%% Velocity triggered windowing

for ii=1:size(Encoder.velTrig,2)
    Waves(ii).p = LFP.xgp(:,:,Encoder.trialTime(ii,3):Encoder.trialTime(ii,4));
    Waves(ii).wt = LFP.wt(:,:,Encoder.trialTime(ii,3):Encoder.trialTime(ii,4));
    %p = arrayfun(@(jj) inpaint_nans(p(:,:,jj)),1:size(p,3));
    Waves(ii).evaluationPoints = find_evaluation_points(Waves(ii).p,pi,0.2);
    %plot_evaluation_points( p, evaluationPoints );
    [Waves(ii).pm,Waves(ii).pd,Waves(ii).dx,Waves(ii).dy] = phase_gradient_complex_multiplication( Waves(ii).p, spacing );
    % Phase gradient directionality 
    [Waves(ii).PGD] = phase_gradient_directionality(Waves(ii).pm,Waves(ii).dx,Waves(ii).dy);
    % Wavelength
    Waves(ii).wl = 1./abs(Waves(ii).pm);
    % Instantaneous speed
    Waves(ii).insts = instantaneous_speed(Waves(ii).wt,Waves(ii).pm);
    Waves(ii).s = speedSpatial(Waves(ii).wt,Waves(ii).pm);
    % Velcity direction unit vector
    [Waves(ii).vx, Waves(ii).vy] = wavefront_direction(Waves(ii).pd,Waves(ii).insts);
    Waves(ii).velDir = atan2(Waves(ii).vy,Waves(ii).vx);
    % divergence calculation
    Waves(ii).source = find_source_points( Waves(ii).evaluationPoints, X, Y, Waves(ii).dx, Waves(ii).dy );
    % phase correlation with distance (\rho_{\phi,d} measure)
    Waves(ii).rho = zeros( 1, length(Waves(ii).evaluationPoints) );
    for jj = 1:length(Waves(ii).evaluationPoints)
        Waves(ii).ph = angle( Waves(ii).p(:,:,Waves(ii).evaluationPoints(jj)) );
        [Waves(ii).rho(jj),~,Waves(ii).D] = phase_correlation_distance( Waves(ii).ph, Waves(ii).source(:,jj), spacing );
    end
    indxBadWave = find(Waves(ii).rho<rhoThres);
    Waves(ii).evaluationPoints(indxBadWave) = [];
    Waves(ii).source(:,indxBadWave) = [];
    Waves(ii).rho(indxBadWave) = [];
    Waves(ii).D(indxBadWave) = [];
    Waves(ii).nWaves = size(Waves(ii).evaluationPoints,2);
    for kk = 1:Waves(ii).nWaves
        st = Waves(ii).evaluationPoints(kk)-2;
        if st<0, st = 0; end
        sp = Waves(ii).evaluationPoints(kk)+2;
        if sp>size(Waves(ii).p,3), sp = Waves(ii).evaluationPoints(kk); end
        Waves(ii).speed(kk) = mean(abs(Waves(ii).insts(:,:,st:sp)),[1 2 3]); % speed in cm/s
        Waves(ii).speed2(kk) = mean(abs(Waves(ii).s(st:sp)),'all'); % speed in cm/s
    end
    Waves(ii).speedpdg = pgdMean(Waves(ii).PGD,Waves(ii).s,0.51);
    Waves(ii).dirpdg = pgdMean(Waves(ii).PGD,Waves(ii).velDir,0.51);
    %plot_vector_field( exp( 1i .* Waves(ii).pd(:,:,100) ), 0 );
    %plot_wave_examples( LFP.xf(:,:,Encoder.trialTime(ii,3):Encoder.trialTime(ii,4)), options, ii, Waves(ii).evaluationPoints, Waves(ii).source, Waves(ii).rho );
end
%%
trialPlot = 5;
plot_wave_examples( LFP.xf(:,:,Encoder.trialTime(trialPlot,3):Encoder.trialTime(trialPlot,4)), ...
    options, trialPlot, Waves(trialPlot).evaluationPoints, Waves(trialPlot).source, Waves(trialPlot).rho,Waves(trialPlot).vx,Waves(trialPlot).vy );
%% Waves accross trials 
speedComb = horzcat(Waves(1:end).speed);
figure();histogram(speedComb,100);
avgSpeed = mean(speedComb);

speedComb2 = horzcat(Waves(1:end).speed2);
figure();histogram(speedComb2,100);
avgSpeed2 = mean(speedComb2);

%% Plotting video of phase 
figure(); hold on;
map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) );axis off;
x = repmat(1:cols,rows,1); % generate x-coordinates
y = repmat(1:rows,cols,1)'; % generate y-coordinates
Ts = 1/1024;
tvel = ((Ts:Ts:300*Ts))*1000;
for i=1:numel(goodTrial.xgp)
    ang = inpaint_nans(angle(goodTrial.xgp(:,:,i)));
%     ang = inpaint_nans(goodTrial.xf(:,:,i));    
    pause(0.1);
    t= num2cell(rad2deg(ang));
    t = cellfun(@num2str, t, 'UniformOutput',false);
%     subplot(2,1,1);
    imagesc(ang); title('Phase map for Electrode Map');
    text(x(:),y(:),t,'HorizontalAlignment','Center');
%     ax2 = axes; set( ax2, 'position', [0.02116    0.80976    0.0484    0.1000] ); axis image
%     [x1,y1] = pol2cart( angle( exp(1i.*linspace(-pi,pi,100)) ), ones( 1, 100 ) );
%     h3 = cline( x1, y1, linspace(-pi,pi,100) ); axis off; set( h3, 'linewidth', 6 )
%     t1 = text( 0, 0, 'GP' );
%     set( t1, 'fontname', 'arial', 'fontsize', 20, 'fontweight', 'bold', 'horizontalalignment', 'center' )
    
%     subplot(2,1,2)
%     plot(Encoder.vel(1+Encoder.velTrig(137)-350:i+Encoder.velTrig(137)-350),'-b','LineWidth',2);xlim([tvel(1) tvel(end)]); ylim([-1 4]);ylabel('Velocity (in cm/s)');xlabel('Time (in ms)');
end

%% Plotting video of LFP amplitude
figure();
x = repmat(1:cols,rows,1); % generate x-coordinates
y = repmat(1:rows,cols,1)'; % generate y-coordinates
tvel = Encoder.timeWindow1;
for i=1:numel(goodTrial.xgp)
    amp = inpaint_nans(goodTrial.xf(:,:,i));    
    pause(0.1);
%     subplot(2,1,1);
    imagesc(amp); title('Phase map for Electrode Map');colorbar;
%     if mod(i,2)==0, continue; end
%     subplot(2,1,2)
%     plot(tvel(1:floor(i/2)),Encoder.vel(1,Encoder.trialTime(trialno,1):Encoder.trialTime(trialno,1)+floor(i/2)),'-b','LineWidth',2)
%     plot(tvel(1:floor(i/2)),'-b','LineWidth',2);
%     xlim([tvel(1) tvel(end)]); ylim([-1 4]);ylabel('Velocity (in cm/s)');xlabel('Time (in s)');
end
%% Finding source points 
% plot_wave_examples( goodTrial.xf, options, jj, evaluationPoints, source, rho );

figure();
for i=1:numel(goodTrial.relTime)
    quiver(dx(:,:,i),dy(:,:,i));
    pause(0.1);
    ylim([0 9]);xlim([0 5]);
end

%% 

% dt = 1 / 1024; T = size(LFP.LFP,2) / Fs; time = dt:dt:T;
time = LFP.times(500:6000);
xw = squeeze(LFP.xf(4,3,500:6000));
xwRaw = squeeze(LFP.LFPdatacube(4,3,500:6000));
xgp1 = squeeze(LFP.xgp(4,3,500:6000));
% main figure
fg1 = figure; hold on; ax1 = gca; 
plot( time, xw, 'linewidth', 2, 'color', 'k' ); h4 = cline( time, xw, [], angle(xgp1) );
xlim([1 2])
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

%% 

waveden = zeros(1,size(LFP.LFP,2));
waveden(evaluationPoints) = 1;
waveden1 = downsample(waveden,2);
figure(),plot(waveden1);


figure();
plot(Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
hold on;
plot(waveden1);



%meanPhase is the mean phase map over a 5ms window starting at the wave
%initiaion time
[rho(jj),~,D,pl] = phase_correlation_distance(meanPhase,source(:,jj), parameters.pixel_spacing);
%fit the phase (p1) vs distance (D) to a line
pfit = polyfit(D,pl,1);
%use phase vs distance to calculate wave speed
k = abs(pfit(1));%rad/mm, WAVE NUMBER, slope of the lin fit (spatial derivative k=dPhase/dx)
IF = nanmean(IF,3);%Hz, INSTANTANEOUS FREQUENCY for each electrode in the small time window. f = dPhase/dt
IF = nanmean(IF(:));%Hz, AVG INSTANTANEOUS FREQUENCY across the grid
w = IF*2*pi;%(1/s)*(rad) = rad/s, AVG INSTANTANEOUS ANGULAR FREQUENCY across the grid. w = 2pi*f
speed(jj) = (w/k).*(1/1000); %(rad/s)/(rad/mm) = (mm/s)*(1m/1000mm) = m/s, speed

