% clear; clc; 
% close all;
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
IntanConcatenate
fpath    = Intan.path; % where on disk do you want the analysis? ideally and SSD...

%% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = 0:Ts:Intan.Tmax-Ts;

%% Removing bad channels from impedance values
[Z,badCh] = readImp(electrode_map,3e6);
% Intan = removeBadCh(Intan,badCh);

%% Loading Encoder Data

[Encoder.pos, Encoder.vel, Encoder.time] = readPos(LFP.times);
figure('Name','Velocity');plot(Encoder.time,Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);

Encoder.vel = Encoder.vel*-1; % Dont run
Encoder.vel(Encoder.vel<0) = 0;
Encoder.vel(end-10000:end) = 0;

%% Find velocity triggers before velovity reaches 2cm/s
Encoder = detectVelTrig(Encoder,2,0.005,100,300);
figure('Name','Velocity');plot(Encoder.time,Encoder.vel,'LineWidth',1.5);ylim([-10 10]);xlabel('Time (in s)');ylabel('Velocity in cm/s');yline([2 -2]);
hold on;xline(Encoder.velTrig/1000);
%% LFP
set(0,'DefaultFigureWindowStyle','normal')
% LFP = fastpreprocess_filtering(flip(Intan.allIntan,1),8192); %Only run for PFF data
LFP = fastpreprocess_filtering(Intan.allIntan,8192);
LFP = bestLFP(LFP);
LFP = bandFilter(LFP,'depth'); % Extract LFPs based on 'depth' or 'single'
% LFPplot(LFP);
LFP = createDataCube(LFP,6,5,Intan.goodChMap); % Creating datacube

%% Generalized Phase 
LFP.xf = bandpass_filter(LFP.LFPdatacube,4,50,4,1024);
LFP.xgp = generalized_phase(LFP.xf,1024,0);
[X,Y] = meshgrid( 1:cols, 1:rows );
%% Velocity triggered windowing
LFP.windowSize = 600;
for ii=1:size(Encoder.velTrig,2)
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
map = colorcet( 'C2' ); colormap( circshift( map, [ 28, 0 ] ) );axis off;
x = repmat(1:5,6,1); % generate x-coordinates
y = repmat(1:6,5,1)'; % generate y-coordinates
for i=1:500
    ang = inpaint_nans(angle(squeeze(LFP.xgp(:,:,i+500))));
    pause(0.1);
    title('n = ' + string(i));
    t= num2cell(rad2deg(ang));
    t = cellfun(@num2str, t, 'UniformOutput',false);
    imagesc(ang);
    text(x(:),y(:),t,'HorizontalAlignment','Center');
end


% dt = 1 / 1024; T = size(LFP.LFP,2) / Fs; time = dt:dt:T;
time = LFP.times(500:6000);
xw = squeeze(xf(4,3,500:6000));
xwRaw = squeeze(LFP.LFPdatacube(4,3,500:6000));
xgp1 = squeeze(xgp(4,3,500:6000));
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
