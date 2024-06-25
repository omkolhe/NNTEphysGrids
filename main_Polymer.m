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
%%
parameters.Fs = 1000;
parameters.ts = 1/parameters.Fs;


load UCLASingle64Ch_chanmap.mat; % load the channel map for the the shank data

finalElectrodeMap = UCLAProbeMap;

IntanConcatenate

% Generating time series from Intan data
Ts = 1/Intan.offsetSample;
Intan.Tmax = Ts * size(Intan.allIntan,2);
Intan.t = 0:Ts:Intan.Tmax-Ts;

%% LFP
set(0,'DefaultFigureWindowStyle','normal')
LFP = fastpreprocess_filtering(Intan.allIntan,5000);

%% Spike analysis 
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
kilosortPath = [pathname,'\kilosort4'];
SpikeClusters = readNPY(fullfile(kilosortPath, 'spike_clusters.npy'));
SpikeSamples = readNPY(fullfile(kilosortPath, 'spike_times.npy'));
SpikeChannel = readNPY(fullfile(kilosortPath,'channel_positions.npy'));
Spikes.SpikeClusters = SpikeClusters; 
Spikes.SpikeSamples = SpikeSamples;
Spikes = clusterSort(Spikes);
Spikes = ISI(Spikes,0.01,20000,0); %Spikes, Interval, Fs

[spikeAmps, spikeDepths, templateDepths, tempAmps, tempsUnW, templateDuration, waveforms, max_site] =...
    spikeTemplatePosition(kilosortPath,ycoords);
for i = 1:length(tempAmps)
    Spikes.Clusters(i).spikeDepth = templateDepths(i);
    Spikes.Clusters(i).channelDepth = max_site(i);
    Spikes.Clusters(i).spikeAmplitude = tempAmps(i);
    Spikes.Clusters(i).waveforms = waveforms(i,:);
    Spikes.Clusters(i).spikeDuration = templateDuration(i)/20000*1000;
end

clear spikeAmps spikeDepths templateDepths tempAmps tempsUnW templateDuration waveforms max_site

%% Seperating spikes for baseline and polymer 
Spikes.ClustersBaseline = Spikes.Clusters;
Spikes.ClustersPolymer = Spikes.Clusters;

for i=1:size(Spikes.Clusters,2)
    indexNotBasline = find(Spikes.Clusters(i).cluster>tbaseline);
    Spikes.ClustersBaseline(i).cluster(indexNotBasline) = [];
    Spikes.ClustersBaseline(i).spikeTime(indexNotBasline) = [];
    Spikes.ClustersPolymer(i).cluster = Spikes.Clusters(i).cluster(indexNotBasline);
    Spikes.ClustersPolymer(i).spikeTime = Spikes.Clusters(i).spikeTime(indexNotBasline);
    if (~isempty(indexNotBasline))
        if (indexNotBasline(1)-1<=0)
            Spikes.ClustersBaseline(i).ISI = [];
            Spikes.ClustersPolymer(i).ISI = Spikes.Clusters(i).ISI;
        else
            Spikes.ClustersBaseline(i).ISI(indexNotBasline-1) = [];
            Spikes.ClustersPolymer(i).ISI = Spikes.Clusters(i).ISI(indexNotBasline-1);
        end
    else
        Spikes.ClustersPolymer(i).ISI = [];
    end
end

%% Get spikes that fire during both baseline and polymer 
Spikes.ClustersOnlyBaseline = [];
Spikes.ClustersOnlyPolymer = [];
Spikes.ClustersBoth = [];
for i=1:size(Spikes.Clusters,2)
   if (~isempty(Spikes.ClustersBaseline(i).cluster) && ~isempty(Spikes.ClustersPolymer(i).cluster))
       Spikes.ClustersBoth = [Spikes.ClustersBoth;i];
   else
       if (isempty(Spikes.ClustersBaseline(i).cluster))
           Spikes.ClustersOnlyPolymer = [Spikes.ClustersOnlyPolymer;i];
       else
           Spikes.ClustersOnlyBaseline = [Spikes.ClustersOnlyBaseline;i];
       end
   end
end

%% Spikes present in both 
Spikes.ClustersBaseline = Spikes.ClustersBaseline(Spikes.ClustersBoth);
Spikes.ClustersPolymer = Spikes.ClustersPolymer(Spikes.ClustersBoth);

%% Plotting changing in firing pattern 

for i = 1:size(Spikes.ClustersBaseline,2)
    Spikes.ISI(i).baseline = mean(Spikes.ClustersBaseline(i).ISI,'all');
    Spikes.ISI(i).polymer = mean(Spikes.ClustersPolymer(i).ISI,'all');
end

ISIBaseline = cell2mat(arrayfun(@(s) s.baseline,Spikes.ISI,'UniformOutput',false));
ISIPolymer = cell2mat(arrayfun(@(s) s.polymer,Spikes.ISI,'UniformOutput',false));

figure,
plotBox2(ISIBaseline,ISIPolymer);
ranksum(ISIBaseline,ISIPolymer)

%% PSD of LFP 
nPoints = 1000;
pxxbaseline = zeros(64,nPoints/2+1);
pxxpolymer = zeros(64,nPoints/2+1);
for i=1:64
    [pxxbaseline(i,:) f] = pwelch(LFP_Baseline.LFP(i,:),1000,300,nPoints,1000);
    [pxxpolymer(i,:) f] = pwelch(LFP_Polymer.LFP(i,:),1000,300,nPoints,1000);
end

figure();
plot(f(1:81),10*log10(pxxbaseline(:,1:81)),'Color', [166/255 14/255 90/255 0.2],'LineWidth',0.5);
hold on;
plot(f(1:81),10*log10(mean(pxxbaseline(:,1:81),1)),'Color', [166/255 14/255 90/255 1],'LineWidth',1.5);
plot(f(1:81),10*log10(pxxpolymer(:,1:81)),'Color', [53/255 189/255 206/255 0.2],'LineWidth',0.5);
plot(f(1:81),10*log10(mean(pxxpolymer(:,1:81),1)),'Color', [53/255 189/255 206/255 1],'LineWidth',1.5);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Average PSD for FA MI Trials');
box off;

figure();
shadedErrorBar(f(1:81),10*log10(pxxbaseline(:,1:81)),{@mean,@std}, 'lineprops', '-r');
hold on;
shadedErrorBar(f(1:81),10*log10(pxxpolymer(:,1:81)),{@mean,@std}, 'lineprops', '-b');
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
box off;