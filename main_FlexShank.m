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
rmpath(genpath('Dependancies/MVGC1'));
%%  PreProcessing
load 32ChFlexShank_chanmap.mat 
parameters.Fs = 1000;
parameters.nElectrodes = 32;
parameters.ts = 1/parameters.Fs;
parameters.windowBeforePull = 1.5; % in seconds
parameters.windowAfterPull = 1.5; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds 
parameters.windowAfterCue = 1.5; % in seconds 
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 1.5; % in seconds 
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated 
parameters.cool = 0; % 1 - cooling , 0 - no cooling 
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.xspacing = 0.03; % Grid spacing in mm between columns 
parameters.yspacing = 0.06; % Grid spacing in mm between rows
parameters.shank = 0; % 1 - if UCLA 64Ch single shank data is present

finalElectrodeMap = electrode_map;

IntanConcatenate
fpath = Intan.path; % where on disk do you want the analysis? ideally and SSD...

% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = 0:Ts:Intan.Tmax-Ts;

%% Removing bad channels from impedance values
[Z,IntanBehaviour.goodChMap,IntanBehaviour.badChMap] = readImp(electrode_map,5e6);
% [Z,IntanBehaviour.goodChMap,IntanBehaviour.badChMap] = readImp(finalElectrodeMap,10e6);
figure('Name','Impedance Test at 1kHz');boxchart(Z); xlabel('n = ' + string(size(Z,1)));ylabel('Impedance (in \Omega)');set(gca,'xticklabel',{[]})
% Intan.badChMap =[21,22];[1,2];[6,31];[5,10,21]; ;2;7];
%Intan = removeBadCh(Intan,Intan.badCh);
IntanBehaviour.badChMap =[21,22];
IntanBehaviour.badChMap =[31,32];
IntanBehaviour.badChMap =[24,25,26];
IntanBehaviour.badChMap =[];
%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(Intan.allIntan,5000);


%% Loading Lever Data 
plotOption = 1;
[Behaviour] = readLever(parameters,LFP.times,plotOption);

%% Reading behaviour data from Intan traces 
plotOption = 1;
IntanBehaviour = readLeverIntan(parameters,LFP.times,Intan.analog_adc_data,Intan.dig_in_data,Behaviour,plotOption);

%% Add trial segmented data to IntanBehaviour Variable
IntanBehaviour = addLFPToBehaviour(IntanBehaviour,LFP,parameters);
% Saving paramters, path, IntanBehaviour to bin file 
savepath = uigetdir(path);
sessionName = [savepath,'/','Day12_BaselineWaves.mat'];
% save(sessionName,"IntanBehaviour","fpath","parameters","-v7.3");
save(sessionName,"IntanBehaviour","fpath","parameters","Waves","-v7.3"); %,"betaWaves","thetaWaves","gammaWaves",

%% PSD of LFP 
nPoints = 4000;
pxx = zeros(parameters.nElectrodes,nPoints/2+1);
for i=1:parameters.nElectrodes
    [pxx(i,:) f] = pwelch(LFP.LFP(i,:),1000,300,nPoints,1000);
end

fmax = 241;
figure();
plot(f(1:fmax),10*log10(pxx(:,1:fmax)),'Color', [166/255 14/255 90/255 0.2],'LineWidth',0.5);
hold on;
plot(f(1:fmax),10*log10(mean(pxx(:,1:fmax),1)),'Color', [166/255 14/255 90/255 1],'LineWidth',1.5);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average PSD for all channels');
box off;

figure();
shadedErrorBar(f(1:fmax),10*log10(pxx(:,1:fmax)),{@mean,@std}, 'lineprops', '-r');
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
box off;
