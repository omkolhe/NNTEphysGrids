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
rows = 8;  % Number of rows of electrodes on the Grid
cols = 4;  % Numevbr of colums of electrodes on the Grid
IntanConcatenate
fpath    = Intan.path; % where on disk do you want the analysis? ideally and SSD...

%% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = Ts:Ts:Intan.Tmax;

%% Removing bad channels from impedance values
[Z,Intan.goodChMap,Intan.badChMap] = readImp(electrode_map,3e6);
figure('Name','Impedance Test at 1kHz');boxchart(Z); xlabel('n = ' + string(size(Z,1)));ylabel('Impedance (in \Omega)');set(gca,'xticklabel',{[]})
% Intan = removeBadCh(Intan,badCh);

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

% Encoder.vel = abs(Encoder.vel);
Encoder.vel(Encoder.vel<0) = 0;

%% Find velocity triggers before velocity reaches 4cm/s
Encoder = detectVelTrig(Encoder,2,0.05,50,100,LFP.times);
figure('Name','Velocity');plot(Encoder.time,Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
hold on;xline(Encoder.velTrig(1,:)/Encoder.fs);
%% Segementing trial windows
windowBeforeTrig = 0.8; % in seconds
windowAfterTrig = 0.2; % in seconds
Encoder = getVelTrigTrials(Encoder,windowBeforeTrig,windowAfterTrig,LFP.Fs);

Encoder.timeWindow1 = -1*windowBeforeTrig:1/Encoder.fs:windowAfterTrig-1/Encoder.fs; % Time series for plotting
Encoder.timeWindow2 = -1*windowBeforeTrig:1/LFP.Fs:windowAfterTrig-1/LFP.Fs; % Time series for plotting
figure('Name','Velocity Triggers trials');
for i=1:1 %Encoder.nTrig
    plot(Encoder.timeWindow1,Encoder.vel(Encoder.trialTime(i,1):Encoder.trialTime(i,2)),'LineWidth',1.5)
    ylim([0 20]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
    hold on;xline(0);
end

%% Generalized Phase 
LFP.xf = bandpass_filter(LFP.LFPdatacube,5,40,4,1000);
LFP.xgp = generalized_phase(LFP.xf,1000,0);
[X,Y] = meshgrid( 1:cols, 1:rows );

%% Wavelet spectrogram
trialno = 2;
goodTrial.xf = LFP.xf(:,:,Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4));
goodTrial.xgp = LFP.xgp(:,:,Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4));
goodTrial.time = LFP.times(Encoder.trialTime(trialno,3):Encoder.trialTime(trialno,4));
goodTrial.relTime = Encoder.timeWindow2;
figure();
for i=1:32
    subplot(rows,cols,i);
    if ismember(i,Intan.badChMap), continue; end
    calCWTSpectogram(squeeze(goodTrial.xf(floor((i-1)/cols)+1,mod(i-1,cols)+1,:)),goodTrial.relTime,1024,20,[1 45],1);
%     [s,f,t,ps,fc,tc] = spectrogram(squeeze(goodTrial.xf(floor((i-1)/cols)+1,mod(i-1,cols)+1,:)),25,16,1024,1024,'yaxis','onesided');
%     imagesc(t,f(10:45),abs(ps(10:45,:)));hold on; xline(301/1024,'-r', 'LineWidth',2);set(gca,'YDir','normal') 
end

% Average spectrogram across all channels
for i=1:rows
    for j=1:cols
        [spectrogramCh((i-1)*cols + j,:,:) ,fwt] = calCWTSpectogram(squeeze(goodTrial.xf(rows,cols,:)),goodTrial.relTime,1024,20,[1 45],0);
    end
end
avgSpectrogramCWT = mean(spectrogramCh,1);
figure();
subplot(2,1,1);
imagesc(goodTrial.relTime,fwt,squeeze(avgSpectrogramCWT(1,:,:)));colormap('jet');set(gca,'YDir','normal');title('Wavelet based Average Spectogram');ylabel('Frequency (Hz)');xlabel('Time (s)');
c=colorbar;ylabel(c, 'Relative Power to white noise','FontSize',10);
subplot(2,1,2);
plot(Encoder.timeWindow1,Encoder.vel(Encoder.trialTime(trialno,1):Encoder.trialTime(trialno,2)),'-r','LineWidth',2);
ylabel('Velocity (cm/s)');xlabel('Time (ms)');

%% PLotting LFP 
figure();
for i=1:32
    subplot(rows,cols,i);
    if ismember(i,Intan.badChMap), continue; end
    plot(goodTrial.relTime,squeeze(goodTrial.xf(floor((i-1)/cols)+1,mod(i-1,cols)+1,:))');xline(0,'-r');
end
%% Velocity triggered windowing

options.subject = 'W'; % this can be 'W' or 'T' (two marmoset subjects)
options.plot = true; % this option turns plots ON or OFF
options.plot_shuffled_examples = false; % example plots w/channels shuffled in space

% LFP.windowSize = 600;
for ii=1:size(Encoder.velTrig,2)
    if Encoder.velTrig(ii) < LFP.windowSize+1 continue; end
    LFP.velTrigTrials(ii,1) = Encoder.velTrig(ii)-LFP.windowSize; % start time 
    LFP.velTrigTrials(ii,2) = Encoder.velTrig(ii)-300; % stop time
    start_time = LFP.velTrigTrials(ii,1)-LFP.windowSize;
    stop_time = LFP.velTrigTrials(ii,2)-300;
    p = LFP.xgp(:,:,LFP.velTrigTrials(ii,1):LFP.velTrigTrials(ii,2));
%     p = arrayfun(@(jj)inpaint_nans(p(:,:,jj)),1:size(p,3));
    evaluationPoints = find_evaluation_points(p,pi,0.2);
%     plot_evaluation_points( p, evaluationPoints );
     [pm,pd,dx,dy] = phase_gradient_complex_multiplication( p, 0.2 );
    
    % divergence calculation
    source = find_source_points( evaluationPoints, X, Y, dx, dy );
    
    % phase correlation with distance (\rho_{\phi,d} measure)
    rho = zeros( 1, length(evaluationPoints) );
    for jj = 1:length(evaluationPoints)
        
        ph = angle( p(:,:,evaluationPoints(jj)) );
        if strcmp( 'W', 'T' ); ph(data(1).mask) = NaN; end
        rho(jj) = phase_correlation_distance( ph, source(:,jj), 0.2 );
        
    end
    % plotting 2 - wave examples
    plot_wave_examples( LFP.xf(:,:,start_time:stop_time), options, ii, evaluationPoints, source, rho );
   
end


%% Plotting video of phase 
figure(); hold on;
% map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) );axis off;
x = repmat(1:cols,rows,1); % generate x-coordinates
y = repmat(1:rows,cols,1)'; % generate y-coordinates
Ts = 1/1024;
tvel = ((Ts:Ts:300*Ts))*1000;
for i=1:numel(goodTrial.xgp)
%     ang = inpaint_nans(angle(goodTrial.xgp(:,:,i)));
    ang = inpaint_nans(goodTrial.xf(:,:,i));    
    pause(0.1);
%     t= num2cell(rad2deg(ang));
%     t = cellfun(@num2str, t, 'UniformOutput',false);
%     subplot(2,1,1);
    imagesc(ang); title('Phase map for Electrode Map');
%     ax2 = axes; set( ax2, 'position', [0.02116    0.80976    0.0484    0.1000] ); axis image
%     [x1,y1] = pol2cart( angle( exp(1i.*linspace(-pi,pi,100)) ), ones( 1, 100 ) );
%     h3 = cline( x1, y1, linspace(-pi,pi,100) ); axis off; set( h3, 'linewidth', 6 )
%     t1 = text( 0, 0, 'GP' );
%     set( t1, 'fontname', 'arial', 'fontsize', 20, 'fontweight', 'bold', 'horizontalalignment', 'center' )
%     text(x(:),y(:),t,'HorizontalAlignment','Center');
%     subplot(2,1,2)
%     plot(Encoder.vel(1+Encoder.velTrig(137)-350:i+Encoder.velTrig(137)-350),'-b','LineWidth',2);xlim([tvel(1) tvel(end)]); ylim([-1 4]);ylabel('Velocity (in cm/s)');xlabel('Time (in ms)');
end

%% Plotting video of LFP amplitude
figure();
x = repmat(1:cols,rows,1); % generate x-coordinates
y = repmat(1:rows,cols,1)'; % generate y-coordinates
tvel = Encoder.timeWindow1;
for i=1:numel(goodTrial.xgp)
%     ang = inpaint_nans(goodTrial.xf(:,:,i));    
    pause(0.1);
%     subplot(2,1,1);
%     imagesc(ang); title('Phase map for Electrode Map');
    if mod(i,2)==0, continue; end
%     subplot(2,1,2)
    plot(tvel(1:floor(i/2)),'-b','LineWidth',2);,Encoder.vel(1,Encoder.trialTime(trialno,1):Encoder.trialTime(trialno,1)+floor(i/2))
%     plot(tvel(1:floor(i/2)),'-b','LineWidth',2);
    %xlim([tvel(1) tvel(end)]); ylim([-1 4]);ylabel('Velocity (in cm/s)');xlabel('Time (in s)');
end
%% Finding source points 

evaluationPoints = find_evaluation_points(goodTrial.xgp,pi,0.2);
plot_evaluation_points(goodTrial.xgp, evaluationPoints);
[X,Y] = meshgrid( 1:cols, 1:rows ); 
[pm,pd,dx,dy] = phase_gradient_complex_multiplication(goodTrial.xgp,0.1);
source = find_source_points(evaluationPoints,X,Y,dx,dy);

rho = zeros( 1, length(evaluationPoints) );
for jj = 1:length(evaluationPoints)
    
    ph = angle( goodTrial.xgp(:,:,evaluationPoints(jj)) );
    rho(jj) = phase_correlation_distance( ph, source(:,jj), 0.1 );
end
plot_wave_examples( goodTrial.xf, options, ii, evaluationPoints, source, rho );

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
